import * as THREE from 'three';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';

export const CELL_URL = '/assets/cell/ayari_cell_cutaway.glb';
export const CELL_GOLD = new THREE.Color(0xc9a96e);

export function collectCellParts(model) {
  const byMarkerKey = new Map();
  const byName = new Map();

  model.traverse((object) => {
    if (!object.userData?.markerKey) return;
    byName.set(object.name, object);
    const markerKey = object.userData.markerKey;
    const group = byMarkerKey.get(markerKey) ?? [];
    group.push(object);
    byMarkerKey.set(markerKey, group);
  });

  return { byMarkerKey, byName };
}

function forEachMaterial(object, callback) {
  object.traverse((child) => {
    if (!child.isMesh) return;
    const materials = Array.isArray(child.material) ? child.material : [child.material];
    materials.forEach((material) => material && callback(material));
  });
}

export function setMarkerHighlight(parts, markerKey, active = true, intensity = 2.2) {
  const targets = parts.byMarkerKey.get(markerKey);
  if (!targets) throw new Error(`Unknown cell markerKey: ${markerKey}`);

  targets.forEach((object) => {
    forEachMaterial(object, (material) => {
      if (!material.userData.ayariOriginalEmissive) {
        material.userData.ayariOriginalEmissive = material.emissive?.clone() ?? new THREE.Color(0x000000);
        material.userData.ayariOriginalEmissiveIntensity = material.emissiveIntensity ?? 1;
      }
      if (!material.emissive) return;
      material.emissive.copy(active ? CELL_GOLD : material.userData.ayariOriginalEmissive);
      material.emissiveIntensity = active ? intensity : material.userData.ayariOriginalEmissiveIntensity;
      material.needsUpdate = true;
    });
  });
}

export function getMarkerBounds(parts, markerKey) {
  const targets = parts.byMarkerKey.get(markerKey);
  if (!targets) throw new Error(`Unknown cell markerKey: ${markerKey}`);
  const bounds = new THREE.Box3();
  targets.forEach((object) => bounds.expandByObject(object, true));
  return bounds;
}

export function focusMarker(camera, controls, parts, markerKey, fitOffset = 1.35) {
  const bounds = getMarkerBounds(parts, markerKey);
  const sphere = bounds.getBoundingSphere(new THREE.Sphere());
  const halfFov = THREE.MathUtils.degToRad(camera.fov * 0.5);
  const distance = (sphere.radius / Math.tan(halfFov)) * fitOffset;
  const direction = camera.position.clone().sub(controls.target).normalize();
  controls.target.copy(sphere.center);
  camera.position.copy(sphere.center).addScaledVector(direction, distance);
  camera.near = Math.max(distance / 100, 0.001);
  camera.far = Math.max(distance * 100, 10);
  camera.updateProjectionMatrix();
  controls.update();
}

export async function loadAyariCell(scene) {
  const gltf = await new GLTFLoader().loadAsync(CELL_URL);
  const model = gltf.scene;
  model.name = 'AYARI_CELL';
  scene.add(model);
  return { model, parts: collectCellParts(model) };
}
