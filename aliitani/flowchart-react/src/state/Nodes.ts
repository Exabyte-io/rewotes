import { proxy } from 'valtio';

let id = 0;
const getId = () => `node_${id++}`;

const initialNodes = [
	{
		id: 'node_start_value',
		type: 'startNode',
		data: {
			label: 'Start',
			startValue: 0,
		},
		position: { x: 250, y: 5 },
	},
	{
		id: 'node_result_value',
		type: 'endNode',
		data: {
			label: 'End of Workflow',
		},
		position: { x: 250, y: 150 },
	},
];

export const Nodes = proxy([...initialNodes]);

export const addNewNode = (node) => {
	Nodes.push(node);
};

export const getNode = (nodeId) => {
	return Nodes.find((node) => node.id === nodeId);
};

export const clearNodes = () => {
	Nodes.splice(0, Nodes.length);
	Nodes.push(...initialNodes);
};
