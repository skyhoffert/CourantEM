/**
 * 3D model view (three.js): the structure being simulated, ports, and after a
 * run a time-domain snapshot of |E| just under the signal layer. Units: mm.
 * Mouse: drag to orbit, wheel to zoom.
 */
import * as THREE from 'three';
import type { Board } from '../import/kicad.ts';
import type { Panel, State } from './store.ts';

const COPPER = 0xd08a3a, SUBSTRATE = 0x3f8f6b, VIA = 0xb8b8b8, PORT = 0x2a7de1;

/** The built-in GCPW as a board model, so one renderer serves both sources. */
export function gcpwAsBoard(len: number): Board {
  const h = 254e-6, w = 500e-6, s = 200e-6, fence = 300e-6;
  const a = w / 2, b = a + s, yf = b + fence, W = yf + 1e-3;
  const rect = (x0: number, y0: number, x1: number, y1: number): [number, number][] => [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
  return {
    layers: [
      { name: 'B.Cu', type: 'copper', thickness: 35e-6, z0: 0, z1: 0 },
      { name: 'core', type: 'dielectric', thickness: h, epsR: 3.66, z0: 0, z1: h },
      { name: 'F.Cu', type: 'copper', thickness: 35e-6, z0: h, z1: h },
    ],
    copperLayers: ['B.Cu', 'F.Cu'],
    tracks: [{ a: [0, 0], b: [len, 0], width: w, layer: 'F.Cu', net: 1 }],
    vias: [],
    polys: [
      { pts: rect(0, -W, len, W), layer: 'B.Cu', net: 2 },
      { pts: rect(0, -W, len, -b), layer: 'F.Cu', net: 2 },
      { pts: rect(0, b, len, W), layer: 'F.Cu', net: 2 },
      // fences drawn as thin copper walls in the side view
    ],
    outline: { min: [0, -W], max: [len, W] },
    nets: new Map([[1, 'RF'], [2, 'GND']]),
    thickness: h,
    walls: [{ y: -yf, z0: 0, z1: h, x0: 0, x1: len }, { y: yf, z0: 0, z1: h, x0: 0, x1: len }],
  } as Board & { walls: { y: number; z0: number; z1: number; x0: number; x1: number }[] };
}

export class ModelPanel implements Panel {
  private el: HTMLElement;
  // three is loaded without type packages in this build; see src/types/vendor.d.ts
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  private renderer?: any;
  private scene = new THREE.Scene();
  private camera = new THREE.PerspectiveCamera(35, 1, 0.01, 1000);
  private root = new THREE.Group();
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  private field?: any;
  private theta = -0.9; private phi = 0.9; private dist = 30;
  private target = new THREE.Vector3();
  private key = '';
  private planes: [number, number] = [0, 0];

  constructor(el: HTMLElement) {
    this.el = el;
    try {
      this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    } catch {
      el.innerHTML = '<p class="muted">3D view needs WebGL.</p>';
      return;
    }
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    el.appendChild(this.renderer.domElement);
    this.scene.add(this.root);
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.75));
    const sun = new THREE.DirectionalLight(0xffffff, 0.9);
    sun.position.set(5, -8, 10);
    this.scene.add(sun);
    // dimmer fill from below so the bottom side isn't lit by ambient alone
    const fill = new THREE.DirectionalLight(0xffffff, 0.5);
    fill.position.set(-5, 8, -10);
    this.scene.add(fill);
    this.bindMouse();
    new ResizeObserver(() => this.render()).observe(el);
  }

  private bindMouse() {
    const c = this.renderer!.domElement;
    let drag: { x: number; y: number } | null = null;
    c.addEventListener('pointerdown', (e: PointerEvent) => { drag = { x: e.clientX, y: e.clientY }; c.setPointerCapture(e.pointerId); });
    c.addEventListener('pointerup', () => { drag = null; });
    c.addEventListener('pointermove', (e: PointerEvent) => {
      if (!drag) return;
      this.theta -= (e.clientX - drag.x) * 0.008;
      // full polar range so the bottom side is reachable; stop short of the poles where up = +z degenerates
      this.phi = Math.min(Math.PI - 0.05, Math.max(0.05, this.phi - (e.clientY - drag.y) * 0.008));
      drag = { x: e.clientX, y: e.clientY };
      this.render();
    });
    c.addEventListener('wheel', (e: WheelEvent) => { e.preventDefault(); this.dist *= Math.exp(e.deltaY * 0.001); this.render(); }, { passive: false });
  }

  update(s: State, prev: State | undefined): void {
    if (!this.renderer) return;
    const board = s.source === 'board' ? s.board?.model : gcpwAsBoard(s.params.len);
    const key = s.source === 'board' ? `b:${s.board?.rev}` : `g:${s.params.len}`;
    this.planes = s.source === 'board' && s.board ? s.board.offsets : [0, 0];
    if (board && key !== this.key) { this.key = key; this.build(board); }
    if (!prev || s.result !== prev.result) this.showField(s);
    this.render();
  }

  private build(b: Board) {
    this.root.clear();
    this.field = undefined;
    const mm = 1e3;
    const [x0, y0] = b.outline.min, [x1, y1] = b.outline.max;
    for (const l of b.layers) {
      if (l.type !== 'dielectric') continue;
      const g = new THREE.BoxGeometry((x1 - x0) * mm, (y1 - y0) * mm, Math.max(1e-3, (l.z1 - l.z0) * mm));
      const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ color: SUBSTRATE, transparent: true, opacity: 0.55 }));
      m.position.set(((x0 + x1) / 2) * mm, ((y0 + y1) / 2) * mm, ((l.z0 + l.z1) / 2) * mm);
      this.root.add(m);
    }
    const cu = new THREE.MeshStandardMaterial({ color: COPPER, metalness: 0.6, roughness: 0.35, side: THREE.DoubleSide });
    const zOf = (layer: string) => (b.layers.find((l) => l.name === layer)?.z0 ?? 0) * mm;
    for (const p of b.polys) {
      const sh = new THREE.Shape(p.pts.map((q) => new THREE.Vector2(q[0] * mm, q[1] * mm)));
      const m = new THREE.Mesh(new THREE.ShapeGeometry(sh), cu);
      m.position.z = zOf(p.layer) + (p.layer === b.copperLayers[0] ? -0.002 : 0.004);
      this.root.add(m);
    }
    for (const t of b.tracks) {
      const dx = (t.b[0] - t.a[0]) * mm, dy = (t.b[1] - t.a[1]) * mm, L = Math.hypot(dx, dy);
      const m = new THREE.Mesh(new THREE.PlaneGeometry(L + t.width * mm, t.width * mm), cu);
      m.position.set(((t.a[0] + t.b[0]) / 2) * mm, ((t.a[1] + t.b[1]) / 2) * mm, zOf(t.layer) + 0.006);
      m.rotation.z = Math.atan2(dy, dx);
      this.root.add(m);
    }
    const viaMat = new THREE.MeshStandardMaterial({ color: VIA, metalness: 0.7, roughness: 0.3 });
    for (const v of b.vias) {
      const h = b.thickness * mm;
      const m = new THREE.Mesh(new THREE.CylinderGeometry(v.drill / 2 * mm, v.drill / 2 * mm, h, 16), viaMat);
      m.rotation.x = Math.PI / 2;
      m.position.set(v.at[0] * mm, v.at[1] * mm, h / 2);
      this.root.add(m);
    }
    const walls = (b as Board & { walls?: { y: number; z0: number; z1: number; x0: number; x1: number }[] }).walls ?? [];
    for (const w of walls) {
      const m = new THREE.Mesh(new THREE.PlaneGeometry((w.x1 - w.x0) * mm, (w.z1 - w.z0) * mm), viaMat);
      m.rotation.x = Math.PI / 2;
      m.position.set(((w.x0 + w.x1) / 2) * mm, w.y * mm, ((w.z0 + w.z1) / 2) * mm);
      this.root.add(m);
    }
    // port planes (at the board edges unless moved into the board)
    const pm = new THREE.MeshBasicMaterial({ color: PORT, transparent: true, opacity: 0.25, side: THREE.DoubleSide });
    for (const x of [x0 + this.planes[0], x1 - this.planes[1]]) {
      const m = new THREE.Mesh(new THREE.PlaneGeometry((y1 - y0) * mm, Math.max(0.5, b.thickness * mm * 3)), pm);
      m.rotation.y = Math.PI / 2; m.rotation.x = Math.PI / 2;
      m.position.set(x * mm, ((y0 + y1) / 2) * mm, b.thickness * mm * 1.5);
      this.root.add(m);
    }
    this.target.set(((x0 + x1) / 2) * mm, ((y0 + y1) / 2) * mm, 0);
    this.dist = 2.3 * Math.max((x1 - x0) * mm, (y1 - y0) * mm, 4);
  }

  private showField(s: State) {
    if (this.field) { this.root.remove(this.field); this.field = undefined; }
    const snap = s.result?.snapshot;
    if (!snap) return;
    const { nx, ny, data } = snap;
    let mx = 0;
    for (let i = 0; i < data.length; i++) mx = Math.max(mx, Math.abs(data[i]));
    const img = new Uint8Array(4 * nx * ny);
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
      const v = Math.min(1, Math.sqrt(Math.abs(data[i * ny + j]) / (mx || 1)));
      const o = 4 * (j * nx + i);
      img[o] = 255 * Math.min(1, 2 * v); img[o + 1] = 255 * Math.max(0, 2 * v - 1); img[o + 2] = 60 * (1 - v); img[o + 3] = 255 * Math.min(1, 1.5 * v);
    }
    const tex = new THREE.DataTexture(img, nx, ny, THREE.RGBAFormat);
    tex.needsUpdate = true;
    const xs = snap.x, ys = snap.y, mm = 1e3;
    const m = new THREE.Mesh(new THREE.PlaneGeometry((xs[xs.length - 1] - xs[0]) * mm, (ys[ys.length - 1] - ys[0]) * mm),
      new THREE.MeshBasicMaterial({ map: tex, transparent: true, side: THREE.DoubleSide, depthWrite: false }));
    m.position.set(((xs[0] + xs[xs.length - 1]) / 2) * mm, ((ys[0] + ys[ys.length - 1]) / 2) * mm, snap.z * mm + 0.01);
    this.field = m;
    this.root.add(m);
  }

  private render() {
    if (!this.renderer) return;
    const w = this.el.clientWidth || 600, h = this.el.clientHeight || 360;
    this.renderer.setSize(w, h, false);
    this.renderer.domElement.style.width = '100%';
    this.renderer.domElement.style.height = `${h}px`;
    this.camera.aspect = w / h;
    this.camera.up.set(0, 0, 1);
    const t = this.target;
    this.camera.position.set(t.x + this.dist * Math.sin(this.phi) * Math.cos(this.theta),
      t.y + this.dist * Math.sin(this.phi) * Math.sin(this.theta), t.z + this.dist * Math.cos(this.phi));
    this.camera.lookAt(t);
    this.camera.updateProjectionMatrix();
    this.renderer.render(this.scene, this.camera);
  }
}
