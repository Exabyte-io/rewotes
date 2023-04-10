import { useCallback } from 'react';
import calculate from '../utils/calculate';

const useUpdateOutputNodes = (nodes, edges, setNodes) => {
    const updateOutputNodes = useCallback(() => {
        setNodes((currentNodes) => {
            // replace all the output nodes in the `nodes` state with new values
            const newNodes = currentNodes.map((node) => {
                if (node.type !== 'outputNode') return node;
                const connectedEdge = edges.find(
                    (edge) => edge.target === node.id
                );
                if (connectedEdge) {
                    const newValue = calculate(
                        currentNodes,
                        edges,
                        connectedEdge.sourceHandle
                    );
                    return { ...node, data: { ...node.data, value: newValue } };
                } else {
                    return node;
                }
            });
            return newNodes;
        });
    }, [nodes, edges, setNodes]);

    return updateOutputNodes;
};
export default useUpdateOutputNodes;