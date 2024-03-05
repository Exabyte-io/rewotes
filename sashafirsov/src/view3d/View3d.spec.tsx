import { render } from '@testing-library/react';
import ReactThreeTestRenderer from '@react-three/test-renderer';

import 'vitest-canvas-mock';

import View3d, { View3dCanvas } from './View3d';
import * as v3mod from './View3d';

import { expect } from 'vitest';

global.ResizeObserver = vi.fn().mockImplementation(() => ({
    observe: vi.fn(),
    unobserve: vi.fn(),
    disconnect: vi.fn()
}));

Object.defineProperty(window, 'matchMedia', {
    writable: true,
    value: vi.fn().mockImplementation(query => ({
        matches: false,
        media: query,
        onchange: null,
        addListener: vi.fn(), // deprecated
        removeListener: vi.fn(), // deprecated
        addEventListener: vi.fn(),
        removeEventListener: vi.fn(),
        dispatchEvent: vi.fn()
    }))
});

describe('View3d', () => {
    afterEach(() => {
        vi.restoreAllMocks()
    })
    it('should render successfully', () => {
        const { baseElement } = render(<View3d />);
        expect(baseElement).toBeTruthy();
    });

    it('should render in canvas', async () => {
        window.localStorage.setItem('xyzdrawing', `{"name":"","slides":[{"comment":"Created by chemcoord http://chemcoord.readthedocs.io/en/latest/","elements":[{"element":"O","x":0,"y":0,"z":0,"sourceLine":3},{"element":"H","x":0.758602,"y":0,"z":0.504284,"sourceLine":4},{"element":"H","x":0.260455,"y":0,"z":-0.872893,"sourceLine":5},{"element":"O","x":3,"y":0.5,"z":0,"sourceLine":6},{"element":"H","x":3.758602,"y":0.5,"z":0.504284,"sourceLine":7},{"element":"H","x":3.260455,"y":0.5,"z":-0.872893,"sourceLine":8}]}]}`);
        window.localStorage.setItem('EditorSelection', `[4,4]`);
        const renderer = await ReactThreeTestRenderer.create(<View3dCanvas />);

        expect(renderer).toBeTruthy();
        expect(renderer.scene.children.length).toEqual(9); // 2 lights, 1 cube, 6 spheres
        expect(renderer.scene.children['8']._fiber.geometry.type).toEqual('SphereGeometry');
        expect(renderer.scene.children.filter(shape => shape._fiber.geometry?.type === 'SphereGeometry').length).toEqual(6);
    });
    it('enlarge on click', async () => {
        window.localStorage.setItem('xyzdrawing', `{"name":"","slides":[{"comment":"Created by chemcoord http://chemcoord.readthedocs.io/en/latest/","elements":[{"element":"O","x":0,"y":0,"z":0,"sourceLine":3},{"element":"H","x":0.758602,"y":0,"z":0.504284,"sourceLine":4},{"element":"H","x":0.260455,"y":0,"z":-0.872893,"sourceLine":5},{"element":"O","x":3,"y":0.5,"z":0,"sourceLine":6},{"element":"H","x":3.758602,"y":0.5,"z":0.504284,"sourceLine":7},{"element":"H","x":3.260455,"y":0.5,"z":-0.872893,"sourceLine":8}]}]}`);
        window.localStorage.setItem('EditorSelection', `[4,4]`);
        const renderer = await ReactThreeTestRenderer.create(<View3dCanvas />);

        expect(renderer).toBeTruthy();
        const prevScale = renderer.scene.children[8].props.scale;
        expect(prevScale).toEqual(1);
        await renderer.fireEvent(renderer.scene.children['8'], 'click');
        const nextScale = renderer.scene.children[8].props.scale;
        expect(nextScale).toEqual(1.5);

    });
});
