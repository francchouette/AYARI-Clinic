import * as THREE from 'three';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';

export const ANATOMY_URL = '/assets/anatomy/ayari_human_anatomy.glb';

export function createGoldMarker(radius = 0.009) {
  const material = new THREE.MeshStandardMaterial({
    color: 0xc9a96e,
    emissive: 0xc9a96e,
    emissiveIntensity: 2.2,
    metalness: 0.28,
    roughness: 0.24,
  });
  const marker = new THREE.Mesh(new THREE.SphereGeometry(radius, 20, 12), material);
  const light = new THREE.PointLight(0xc9a96e, 0.45, 0.14);
  marker.add(light);
  return marker;
}

export function collectAnatomyLocators(model) {
  const locators = new Map();
  model.traverse((object) => {
    if (object.userData?.nodeType === 'locator') {
      locators.set(object.userData.markerKey ?? object.name, object);
    }
  });
  return locators;
}

export function attachMarker(locators, markerKey, marker = createGoldMarker()) {
  const locator = locators.get(markerKey);
  if (!locator) {
    throw new Error(`Unknown anatomy marker: ${markerKey}`);
  }
  locator.add(marker);
  return marker;
}

export async function loadAyariAnatomy(scene) {
  const loader = new GLTFLoader();
  const gltf = await loader.loadAsync(ANATOMY_URL);
  const model = gltf.scene;
  model.name = 'AYARI_ANATOMY';
  scene.add(model);

  return {
    model,
    locators: collectAnatomyLocators(model),
    attachMarker(markerKey, marker) {
      return attachMarker(this.locators, markerKey, marker);
    },
  };
}
