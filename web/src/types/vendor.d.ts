// Minimal module shims so the app type-checks without @types packages installed.
// With `npm install`, uplot ships its own types and @types/three can replace the shim.
declare module 'three';
declare module 'uplot';
declare module '*.css' {
  const text: string;
  export default text;
}

/** CalVer of this build (esbuild define from /VERSION). */
declare const __COURANT_VERSION__: string;
