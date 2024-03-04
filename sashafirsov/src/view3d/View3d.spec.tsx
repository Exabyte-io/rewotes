import { render } from '@testing-library/react';
import 'vitest-canvas-mock'

import View3d from './View3d';
import App from '../app/app';

global.ResizeObserver = vi.fn().mockImplementation(() => ({
    observe: vi.fn(),
    unobserve: vi.fn(),
    disconnect: vi.fn(),
}))

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
        dispatchEvent: vi.fn(),
    })),
})

describe('App', () => {
    it('should render successfully', () => {
        const { baseElement } = render(<View3d />);
        expect(baseElement).toBeTruthy();
    });

    it('should have a tooltip', () => {
        const { queryByTestId } = render(<View3d />);
        expect(queryByTestId('view3d-tooltip')).toBeTruthy();
    });

    it('should have a tooltip', () => {
        const { queryByTestId } = render(<View3d />);
        const canvas = queryByTestId('view3d-canvas');
        expect(canvas).toBeTruthy();
    });
});
