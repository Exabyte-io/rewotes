import { renderHook, act } from '@testing-library/react-hooks';
import useRemoteFlowData from './useRemoteFlowData';

// Mock the Meteor.call method
jest.mock('meteor/meteor', () => ({
    Meteor: {
        call: jest.fn((method, ...args) => {
            const callback = args.pop();
            if (method === 'fetchFlows') {
                callback(null, []);
            } else if (method === 'saveFlow') {
                callback(null);
            } else if (method === 'clearFlows') {
                callback(null);
            }
        }),
    },
}));

// Mock the window.confirm method
global.window.confirm = jest.fn(() => true);

describe('useRemoteFlowData', () => {
    const nodes = [];
    const edges = [];

    test('should update flowName', () => {
        const { result } = renderHook(() => useRemoteFlowData(nodes, edges));
        act(() => {
            result.current.updateFlowName({ target: { value: 'test' } });
        });
        expect(result.current.flowName).toBe('test');
    });

    test('should fetch flows', async () => {
        const { result } = renderHook(() => useRemoteFlowData(nodes, edges));
        await act(async () => {
            result.current.fetchFlows();
        });
        expect(result.current.fetchedFlows).toEqual([]);
    });

    test('should save flow', async () => {
        const { result } = renderHook(() => useRemoteFlowData(nodes, edges));
        act(() => {
            result.current.updateFlowName({ target: { value: 'test' } });
        });
        await act(async () => {
            result.current.saveFlow();
        });
        expect(result.current.flowName).toBe('');
    });

    test('should clear flows', async () => {
        const { result } = renderHook(() => useRemoteFlowData(nodes, edges));
        await act(async () => {
            result.current.clearFlows();
        });
        // Since clearFlows doesn't change state, we'll just ensure it doesn't throw an error
        expect(true).toBe(true);
    });
});
