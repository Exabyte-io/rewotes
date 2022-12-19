import { proxy } from 'valtio';

let id = 0;

const getId = () => `node_${id++}`;

const initialNodes = {
	node_start_value: {
		id: 'node_start_value',
		type: 'startNode',
		data: {
			label: 'Start',
			startValue: 0,
		},
		position: { x: 250, y: 5 },
	},
	node_result_value: {
		id: 'node_result_value',
		type: 'endNode',
		data: {
			label: 'End of Workflow',
			outputValue: 0,
		},
		position: { x: 250, y: 150 },
	},
};

export const Nodes = proxy(initialNodes);

export const addNewNode = (type, position, data) => {
	const nodeId = getId();
	Nodes[nodeId] = {
		id: nodeId,
		type: type,
		position: position,
		data: data,
	};
};

export const updateNodePosition = (position, nodeId) => {
	if (!position) {
		return;
	}
	Nodes[nodeId].position = position;
};

export const getNode = (nodeId) => {
	return Nodes[nodeId];
};

export const clearWorkflow = () => {
	window.location.reload();
};

const initialJsonData = {
	startValue: 0,
	outputValue: 0,
	workflow: [],
	isStartNodeConnected: false,
	isWorkflowValid: false,
	isExecutionDone: false,
};

export const JsonDataState = proxy(initialJsonData);

export const updateObjectValueInJsonViewer = (value, obj) => {
	JsonDataState[obj] = value;
};

export const addNewOperation = (source, target, operation, v) => {
	JsonDataState.workflow.push({ source, target, operation, v });
};

export const updateOperationByNode = (node) => {
	JsonDataState.workflow = JsonDataState.workflow.map((op) => {
		if (op.target === node.id) {
			op.v = node.data.value;
		}
		return op;
	});
};

export const updateConditionByNode = (node, property) => {
	JsonDataState.workflow = JsonDataState.workflow.map((c) => {
		if (c.target === node.id) {
			c[property] = node.data[property];
		}
		return c;
	});
};

export const updateTruthStateConditionByNode = (
	nodeId,
	sourceHandle,
	value
) => {
	JsonDataState.workflow = JsonDataState.workflow.map((c) => {
		if (c.target === nodeId) {
			c[sourceHandle] = value;
		}
		return c;
	});
};

export const updateConnectionStateForConditionalNodes = (p, e) => {
	const nodeExists = JsonDataState.workflow.find((w) => p.source === w.target);
	if (nodeExists) {
		JsonDataState.workflow = JsonDataState.workflow.map((w) => {
			if (nodeExists.target === w.target) {
				return {
					...w,
					falseState: nodeExists.falseState
						? nodeExists.falseState
						: 'node_result_value',
					truthState: nodeExists.truthState
						? nodeExists.truthState
						: 'node_result_value',
				};
			}
			return w;
		});
	} else {
		const target = p.source;
		const valueToGoInTruthFalseState = p.target;
		const location = p.sourceHandle;

		JsonDataState.workflow.push({
			target: target,
			condition: '<',
			value: 0,
			operation: 'conditionals',
			[location]: valueToGoInTruthFalseState,
		});
	}
};

export const addNewCondition = (source, target, condition, value, type) => {
	const nodeExists = JsonDataState.workflow.find((w) => w.target === target);

	if (nodeExists) {
		JsonDataState.workflow = JsonDataState.workflow.map((w) => {
			if (nodeExists.target === w.target) {
				return {
					...w,
					source: source,
					condition: condition,
					type: type,
					value: value,
					falseState: nodeExists.falseState
						? nodeExists.falseState
						: 'node_result_value',
					truthState: nodeExists.truthState
						? nodeExists.truthState
						: 'node_result_value',
				};
			}
			return w;
		});
	} else {
		JsonDataState.workflow.push({
			source,
			target,
			condition: condition ? condition : '<',
			value: value ? value : 0,
			operation: type,
			truthState: 'node_result_value',
			falseState: 'node_result_value',
		});
	}
};
