/**
 * WGSL kernels for the f32 WebGPU backend. Index layout, update ranges and
 * operation order mirror solver/src/sim.c exactly; the C core builds the plan
 * (coefficients, CPML profiles, sources, probes) and this code only time-steps.
 *
 * Storage: E, H, CA, CB are component-major arrays of 3N f32 (comp * N + q).
 * SP packs the per-axis 1D arrays; its offsets are in the uniform block.
 */
export const fdtdWGSL = /* wgsl */ `
struct Params {
  nx: u32, ny: u32, nz: u32, N: u32,
  si: u32, sj: u32, dth: f32, pad0: u32,
  // offsets into SP: ide (n+1), idh (n), be, ce (n+1), bh, ch (n), idd_raw (n+1), id_raw (n) per axis
  ide: vec4<u32>, idh: vec4<u32>, be: vec4<u32>, ce: vec4<u32>,
  bh: vec4<u32>, ch: vec4<u32>, iddr: vec4<u32>, idr: vec4<u32>,
  // CPML slabs per axis: (lo start, lo count, hi start, hi count) in nodes for E, cells for H
  slabE: array<vec4<u32>, 3>,
  slabH: array<vec4<u32>, 3>,
  // offsets into PSI of psiE[c][w] (slot 2c+w) and psiH[c][w] (slot 6+2c+w); compact
  // layout: one entry per thread t of the matching CPML kernel (plane-major, as C)
  psiOff: array<vec4<u32>, 3>,
};

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read_write> E: array<f32>;
@group(0) @binding(2) var<storage, read_write> H: array<f32>;
@group(0) @binding(3) var<storage, read> CA: array<f32>;
@group(0) @binding(4) var<storage, read> CB: array<f32>;
@group(0) @binding(5) var<storage, read> SP: array<f32>;
@group(0) @binding(6) var<storage, read_write> PSI: array<f32>;   // compact CPML auxiliaries, see psiOff

fn ijk(q: u32) -> vec3<u32> {
  let i = q / P.si;
  let r = q - i * P.si;
  let j = r / P.sj;
  return vec3<u32>(i, j, r - j * P.sj);
}

@compute @workgroup_size(256)
fn updateH(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let q = gid.x + gid.y * nwg.x * 256u;
  if (q >= P.N) { return; }
  let c = ijk(q);
  let N = P.N; let si = P.si; let sj = P.sj;
  // Hx: i <= nx, j < ny, k < nz
  if (c.y < P.ny && c.z < P.nz) {
    let iy = SP[P.idh.y + c.y]; let iz = SP[P.idh.z + c.z];
    H[q] -= P.dth * ((E[2u * N + q + sj] - E[2u * N + q]) * iy - (E[N + q + 1u] - E[N + q]) * iz);
  }
  // Hy: i < nx, k < nz
  if (c.x < P.nx && c.z < P.nz) {
    let iz = SP[P.idh.z + c.z]; let ix = SP[P.idh.x + c.x];
    H[N + q] -= P.dth * ((E[q + 1u] - E[q]) * iz - (E[2u * N + q + si] - E[2u * N + q]) * ix);
  }
  // Hz: i < nx, j < ny
  if (c.x < P.nx && c.y < P.ny) {
    let ix = SP[P.idh.x + c.x]; let iy = SP[P.idh.y + c.y];
    H[2u * N + q] -= P.dth * ((E[N + q + si] - E[N + q]) * ix - (E[q + sj] - E[q]) * iy);
  }
}

@compute @workgroup_size(256)
fn updateE(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let q = gid.x + gid.y * nwg.x * 256u;
  if (q >= P.N) { return; }
  let c = ijk(q);
  let N = P.N; let si = P.si; let sj = P.sj;
  // Ex: i < nx, 1 <= j < ny, 1 <= k < nz
  if (c.x < P.nx && c.y >= 1u && c.y < P.ny && c.z >= 1u && c.z < P.nz) {
    let iy = SP[P.ide.y + c.y]; let iz = SP[P.ide.z + c.z];
    E[q] = CA[q] * E[q] + CB[q] * ((H[2u * N + q] - H[2u * N + q - sj]) * iy - (H[N + q] - H[N + q - 1u]) * iz);
  }
  // Ey: 1 <= i < nx, j < ny, 1 <= k < nz
  if (c.x >= 1u && c.x < P.nx && c.y < P.ny && c.z >= 1u && c.z < P.nz) {
    let iz = SP[P.ide.z + c.z]; let ix = SP[P.ide.x + c.x];
    E[N + q] = CA[N + q] * E[N + q] + CB[N + q] * ((H[q] - H[q - 1u]) * iz - (H[2u * N + q] - H[2u * N + q - si]) * ix);
  }
  // Ez: 1 <= i < nx, 1 <= j < ny, k < nz
  if (c.x >= 1u && c.x < P.nx && c.y >= 1u && c.y < P.ny && c.z < P.nz) {
    let ix = SP[P.ide.x + c.x]; let iy = SP[P.ide.y + c.y];
    E[2u * N + q] = CA[2u * N + q] * E[2u * N + q] + CB[2u * N + q] * ((H[N + q] - H[N + q - si]) * ix - (H[q] - H[q - sj]) * iy);
  }
}

// ---- CPML -----------------------------------------------------------------
// Thread t covers one (slab index, other two indices) triple; the slab index is
// mapped onto the low or high layer of the given axis.
fn slabIndex(s: vec4<u32>, t: u32) -> u32 {
  if (t < s.y) { return s.x + t; }
  return s.z + (t - s.y);
}

fn psiSlot(n: u32) -> u32 { return P.psiOff[n / 4u][n % 4u]; }
fn psiE(c: u32, w: u32, t: u32) -> u32 { return psiSlot(2u * c + w) + t; }
fn psiH(c: u32, w: u32, t: u32) -> u32 { return psiSlot(6u + 2u * c + w) + t; }

@compute @workgroup_size(256)
fn pmlHx(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let t = gid.x + gid.y * nwg.x * 256u;
  let s = P.slabH[0];
  let per = (P.ny + 1u) * (P.nz + 1u);
  if (t >= (s.y + s.w) * per) { return; }
  let i = slabIndex(s, t / per);
  let r = t % per; let j = r / (P.nz + 1u); let k = r % (P.nz + 1u);
  let q = i * P.si + j * P.sj + k;
  let N = P.N;
  let b = SP[P.bh.x + i]; let cc = SP[P.ch.x + i]; let id = SP[P.idr.x + i];
  if (k < P.nz) {
    let p = psiH(1u, 1u, t);
    PSI[p] = b * PSI[p] + cc * (E[2u * N + q + P.si] - E[2u * N + q]) * id;
    H[N + q] += P.dth * PSI[p];
  }
  if (j < P.ny) {
    let p = psiH(2u, 0u, t);
    PSI[p] = b * PSI[p] + cc * (E[N + q + P.si] - E[N + q]) * id;
    H[2u * N + q] -= P.dth * PSI[p];
  }
}

@compute @workgroup_size(256)
fn pmlHy(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let t = gid.x + gid.y * nwg.x * 256u;
  let s = P.slabH[1];
  let per = (P.nx + 1u) * (P.nz + 1u);
  if (t >= (s.y + s.w) * per) { return; }
  let j = slabIndex(s, t / per);
  let r = t % per; let i = r / (P.nz + 1u); let k = r % (P.nz + 1u);
  let q = i * P.si + j * P.sj + k;
  let N = P.N;
  let b = SP[P.bh.y + j]; let cc = SP[P.ch.y + j]; let id = SP[P.idr.y + j];
  if (i < P.nx) {
    let p = psiH(2u, 1u, t);
    PSI[p] = b * PSI[p] + cc * (E[q + P.sj] - E[q]) * id;
    H[2u * N + q] += P.dth * PSI[p];
  }
  if (k < P.nz) {
    let p = psiH(0u, 0u, t);
    PSI[p] = b * PSI[p] + cc * (E[2u * N + q + P.sj] - E[2u * N + q]) * id;
    H[q] -= P.dth * PSI[p];
  }
}

@compute @workgroup_size(256)
fn pmlHz(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let t = gid.x + gid.y * nwg.x * 256u;
  let s = P.slabH[2];
  let per = (P.nx + 1u) * (P.ny + 1u);
  if (t >= (s.y + s.w) * per) { return; }
  let k = slabIndex(s, t / per);
  let r = t % per; let i = r / (P.ny + 1u); let j = r % (P.ny + 1u);
  let q = i * P.si + j * P.sj + k;
  let N = P.N;
  let b = SP[P.bh.z + k]; let cc = SP[P.ch.z + k]; let id = SP[P.idr.z + k];
  if (j < P.ny) {
    let p = psiH(0u, 1u, t);
    PSI[p] = b * PSI[p] + cc * (E[N + q + 1u] - E[N + q]) * id;
    H[q] += P.dth * PSI[p];
  }
  if (i < P.nx) {
    let p = psiH(1u, 0u, t);
    PSI[p] = b * PSI[p] + cc * (E[q + 1u] - E[q]) * id;
    H[N + q] -= P.dth * PSI[p];
  }
}

@compute @workgroup_size(256)
fn pmlEx(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let t = gid.x + gid.y * nwg.x * 256u;
  let s = P.slabE[0];
  let per = (P.ny + 1u) * (P.nz + 1u);
  if (t >= (s.y + s.w) * per) { return; }
  let i = slabIndex(s, t / per);
  let r = t % per; let j = r / (P.nz + 1u); let k = r % (P.nz + 1u);
  let q = i * P.si + j * P.sj + k;
  let N = P.N;
  let b = SP[P.be.x + i]; let cc = SP[P.ce.x + i]; let id = SP[P.iddr.x + i];
  if (j < P.ny && k >= 1u && k < P.nz) {
    let p = psiE(1u, 1u, t);
    PSI[p] = b * PSI[p] + cc * (H[2u * N + q] - H[2u * N + q - P.si]) * id;
    E[N + q] -= CB[N + q] * PSI[p];
  }
  if (j >= 1u && j < P.ny && k < P.nz) {
    let p = psiE(2u, 0u, t);
    PSI[p] = b * PSI[p] + cc * (H[N + q] - H[N + q - P.si]) * id;
    E[2u * N + q] += CB[2u * N + q] * PSI[p];
  }
}

@compute @workgroup_size(256)
fn pmlEy(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let t = gid.x + gid.y * nwg.x * 256u;
  let s = P.slabE[1];
  let per = (P.nx + 1u) * (P.nz + 1u);
  if (t >= (s.y + s.w) * per) { return; }
  let j = slabIndex(s, t / per);
  let r = t % per; let i = r / (P.nz + 1u); let k = r % (P.nz + 1u);
  let q = i * P.si + j * P.sj + k;
  let N = P.N;
  let b = SP[P.be.y + j]; let cc = SP[P.ce.y + j]; let id = SP[P.iddr.y + j];
  if (i >= 1u && i < P.nx && k < P.nz) {
    let p = psiE(2u, 1u, t);
    PSI[p] = b * PSI[p] + cc * (H[q] - H[q - P.sj]) * id;
    E[2u * N + q] -= CB[2u * N + q] * PSI[p];
  }
  if (i < P.nx && k >= 1u && k < P.nz) {
    let p = psiE(0u, 0u, t);
    PSI[p] = b * PSI[p] + cc * (H[2u * N + q] - H[2u * N + q - P.sj]) * id;
    E[q] += CB[q] * PSI[p];
  }
}

@compute @workgroup_size(256)
fn pmlEz(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let t = gid.x + gid.y * nwg.x * 256u;
  let s = P.slabE[2];
  let per = (P.nx + 1u) * (P.ny + 1u);
  if (t >= (s.y + s.w) * per) { return; }
  let k = slabIndex(s, t / per);
  let r = t % per; let i = r / (P.ny + 1u); let j = r % (P.ny + 1u);
  let q = i * P.si + j * P.sj + k;
  let N = P.N;
  let b = SP[P.be.z + k]; let cc = SP[P.ce.z + k]; let id = SP[P.iddr.z + k];
  if (i < P.nx && j >= 1u && j < P.ny) {
    let p = psiE(0u, 1u, t);
    PSI[p] = b * PSI[p] + cc * (H[N + q] - H[N + q - 1u]) * id;
    E[q] -= CB[q] * PSI[p];
  }
  if (i >= 1u && i < P.nx && j < P.ny) {
    let p = psiE(1u, 0u, t);
    PSI[p] = b * PSI[p] + cc * (H[q] - H[q - 1u]) * id;
    E[N + q] += CB[N + q] * PSI[p];
  }
}
`;

/** Sources and probes: a second module so the main one stays within 8 storage bindings. */
export const ioWGSL = /* wgsl */ `
struct IO {
  nsrc: u32, nprobe: u32, p0: u32, p1: u32,   // probe range [p0, p1) for this dispatch
};
@group(0) @binding(0) var<uniform> U: IO;
@group(0) @binding(1) var<storage, read_write> E: array<f32>;
@group(0) @binding(2) var<storage, read> H: array<f32>;
@group(0) @binding(3) var<storage, read> SRC: array<vec4<u32>>;     // (index, coef bits, type, 0)
@group(0) @binding(4) var<storage, read> WAVE: array<vec2<f32>>;    // per step: (lumped, mode)
@group(0) @binding(5) var<storage, read> PROBE: array<vec4<u32>>;   // per probe: (start, count, field, 0)
@group(0) @binding(6) var<storage, read> PENT: array<vec2<u32>>;    // per entry: (index, weight bits)
@group(0) @binding(7) var<storage, read_write> STATE: array<u32>;   // [0] = local step; then series (f32 bits)

@compute @workgroup_size(64)
fn inject(@builtin(global_invocation_id) gid: vec3<u32>) {
  let e = gid.x;
  if (e >= U.nsrc) { return; }
  let s = SRC[e];
  let w = WAVE[STATE[0]];
  let v = select(w.x, w.y, s.z == 1u);
  E[s.x] += bitcast<f32>(s.y) * v;
}

@compute @workgroup_size(64)
fn probes(@builtin(global_invocation_id) gid: vec3<u32>) {
  let p = U.p0 + gid.x;
  if (p >= U.p1) { return; }
  let d = PROBE[p];
  var acc = 0.0;
  for (var e = 0u; e < d.y; e++) {
    let en = PENT[d.x + e];
    let v = select(E[en.x], H[en.x], d.z == 1u);
    acc += v * bitcast<f32>(en.y);
  }
  STATE[1u + STATE[0] * U.nprobe + p] = bitcast<u32>(acc);
}

@compute @workgroup_size(1)
fn tick() { STATE[0] = STATE[0] + 1u; }
`;

/** Debye ADE (see solver/src/dispersive.c). One thread per dispersive edge. */
export const dispWGSL = /* wgsl */ `
struct DU { nd: u32, npole: u32, invdt: f32, pad: u32, k: array<vec4<f32>, 4> };
@group(0) @binding(0) var<uniform> U: DU;
@group(0) @binding(1) var<storage, read_write> E: array<f32>;
@group(0) @binding(2) var<storage, read> DL: array<vec2<u32>>;      // (global E index, cb bits)
@group(0) @binding(3) var<storage, read> DB: array<f32>;            // beta, pole-major
@group(0) @binding(4) var<storage, read_write> DJ: array<f32>;      // J, pole-major
@group(0) @binding(5) var<storage, read_write> DO: array<f32>;      // E before the update

fn kp(p: u32) -> f32 { return U.k[p / 4u][p % 4u]; }

@compute @workgroup_size(256)
fn dispSave(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let e = gid.x + gid.y * nwg.x * 256u;
  if (e >= U.nd) { return; }
  DO[e] = E[DL[e].x];
}

@compute @workgroup_size(256)
fn dispCorrect(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let e = gid.x + gid.y * nwg.x * 256u;
  if (e >= U.nd) { return; }
  var acc = 0.0;
  for (var p = 0u; p < U.npole; p++) { acc += 0.5 * (1.0 + kp(p)) * DJ[p * U.nd + e]; }
  let d = DL[e];
  E[d.x] -= bitcast<f32>(d.y) * acc;
}

@compute @workgroup_size(256)
fn dispUpdate(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let e = gid.x + gid.y * nwg.x * 256u;
  if (e >= U.nd) { return; }
  let dE = (E[DL[e].x] - DO[e]) * U.invdt;
  for (var p = 0u; p < U.npole; p++) {
    let q = p * U.nd + e;
    DJ[q] = kp(p) * DJ[q] + DB[q] * dE;
  }
}
`;

/** Surface-impedance copper (see solver/src/sibc.c). One thread per surface edge. */
export const sibcWGSL = /* wgsl */ `
struct SU { ns: u32, nterm: u32, pad0: u32, pad1: u32, c: array<vec4<f32>, 4>, e: array<vec4<f32>, 4> };
@group(0) @binding(0) var<uniform> U: SU;
@group(0) @binding(1) var<storage, read_write> E: array<f32>;
@group(0) @binding(2) var<storage, read> H: array<f32>;
@group(0) @binding(3) var<storage, read> SL: array<u32>;   // per edge 12 words: idx, 4 H idx, 4 weight bits, pad
@group(0) @binding(4) var<storage, read_write> SY: array<f32>;

@compute @workgroup_size(256)
fn sibc(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
  let s = gid.x + gid.y * nwg.x * 256u;
  if (s >= U.ns) { return; }
  let b = 12u * s;
  var J = 0.0;
  for (var f = 0u; f < 4u; f++) { J += bitcast<f32>(SL[b + 5u + f]) * H[SL[b + 1u + f]]; }
  var Ev = 0.0;
  for (var k = 0u; k < U.nterm; k++) {
    let ek = U.e[k / 4u][k % 4u];
    let q = s * U.nterm + k;
    SY[q] = ek * SY[q] + (1.0 - ek) * J;
    Ev += U.c[k / 4u][k % 4u] * (J - SY[q]);
  }
  E[SL[b]] = Ev;
}
`;
