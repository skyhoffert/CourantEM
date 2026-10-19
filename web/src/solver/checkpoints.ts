/**
 * Run checkpoints in IndexedDB (works in workers). A checkpoint is the C core's
 * cem_checkpoint_save() buffer, keyed by a hash of the run request, so an
 * interrupted run (closed tab, crash, reload) resumes when the same run is
 * started again. Deleted when the run finishes.
 */
const DB = 'courant-em', STORE = 'checkpoints';

function open(): Promise<IDBDatabase> {
  return new Promise((res, rej) => {
    const r = indexedDB.open(DB, 1);
    r.onupgradeneeded = () => r.result.createObjectStore(STORE);
    r.onsuccess = () => res(r.result);
    r.onerror = () => rej(r.error);
  });
}

async function tx<T>(mode: IDBTransactionMode, fn: (s: IDBObjectStore) => IDBRequest<T>): Promise<T> {
  const db = await open();
  // resolve on transaction completion (committed), not on request success
  return new Promise((res, rej) => {
    const t = db.transaction(STORE, mode, { durability: 'strict' } as IDBTransactionOptions);
    const r = fn(t.objectStore(STORE));
    t.oncomplete = () => { db.close(); res(r.result); };
    t.onerror = t.onabort = () => { db.close(); rej(t.error ?? r.error); };
  });
}

export async function runKey(request: unknown): Promise<string> {
  // the backend and thread count don't change the physics: a checkpoint resumes on any of them
  const text = JSON.stringify(request, (k, v) => (k === 'wasmUrl' || k === 'threads' || k === 'backend' ? undefined : v));
  const h = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text));
  return Array.from(new Uint8Array(h), (b) => b.toString(16).padStart(2, '0')).join('');
}

/** Store a checkpoint; resolves false (never throws) if storage is unavailable or full. */
export const saveCheckpoint = (key: string, buf: Float64Array) =>
  tx('readwrite', (s) => s.put(buf, key)).then(() => true, (e) => { console.warn('checkpoint save failed', e); return false; });
export const loadCheckpoint = (key: string) =>
  tx<Float64Array | undefined>('readonly', (s) => s.get(key)).catch((e) => { console.warn('checkpoint load failed', e); return undefined; });
export const dropCheckpoint = (key: string) => tx('readwrite', (s) => s.delete(key)).catch(() => undefined);
