import { useSnapshot } from 'valtio';
import { JsonDataState, clearWorkflow } from '../../state/AppState';

const ARITHMETIC_OPERATIONS = ['increment', 'decrement', 'multiply', 'divide'];

/**
 * This helper function handles the arithmetic operations that run against the output value.
 */
const handleArithmeticOperation = (node, currentValue) => {
	let result = currentValue;
	switch (node.operation) {
		case 'increment':
			result += node.v;
			break;
		case 'decrement':
			result -= node.v;
			break;
		case 'multiply':
			result *= node.v;
			break;
		case 'divide':
			result /= node.v;
			break;
		default:
			break;
	}
	return result;
};

/***
 * This helper function handles the conditional logic that runs against the running output value.
 * It runs the conditional logic based on the conditional operators compared against a value;
 */
const handleCondition = (node, currentValue) => {
	let conditionResult = false;
	const condition = node.condition;

	switch (condition) {
		case '<':
			conditionResult = currentValue < node.value;
			break;
		case '>':
			conditionResult = currentValue > node.value;
			break;
		case '<=':
			conditionResult = currentValue <= node.value;
			break;
		case '>=':
			conditionResult = currentValue >= node.value;
			break;
		case '==':
			conditionResult = currentValue == node.value;
			break;
		case '!=':
			conditionResult = currentValue != node.value;
			break;
		default:
			break;
	}

	return conditionResult;
};

const getCurrentNode = (workflow, sourceNode, targetNode) => {
	sourceNode = sourceNode === null ? 'node_start_value' : sourceNode;
	if (targetNode !== null) {
		return workflow.find(
			(w) => w.source === sourceNode && w.target === targetNode
		);
	}

	return workflow.find((w) => {
		return w.source === sourceNode;
	});
};

/***
 * This function runs the entire workflow that's executed. It pulls the workflow data from the AppState.
 * It utilizes two helper functions to handle arithmetic operations and conditional logic.
 * Once the execution is done, it returns a result back that is then updated to
 * the AppState and displayed on the JSON Viewer panel.
 */
const runWorkflow = (json) => {
	const workflow = Array.from(json.workflow);
	let result = json.startValue;

	let nodeRunner = null;
	let nodeTarget = null;

	while (nodeRunner !== 'node_result_value') {
		const currentNode = getCurrentNode(workflow, nodeRunner, nodeTarget);

		if (ARITHMETIC_OPERATIONS.includes(currentNode.operation)) {
			result = handleArithmeticOperation(currentNode, result);
			nodeRunner = currentNode.target;
			nodeTarget = null;
		} else if (currentNode.operation === 'conditionals') {
			const condition = handleCondition(currentNode, result);
			nodeRunner = currentNode.target;
			nodeTarget = condition ? currentNode.truthState : currentNode.falseState;
		} else if (currentNode.operation === 'end_workflow') {
			nodeRunner = 'node_result_value';
			nodeTarget = null;
		}
	}

	return result;
};
/***
 * This component holds the workflow actions that gets triggered against the flowchart.
 * - ExecuteWorkflow - executes the workflow based on the nodes(elements) rendered into the flowchart.
 * - ResetWorkflow - resets the workflow and clears out any state.
 */
const WorkflowActions = () => {
	const json = useSnapshot(JsonDataState);

	const onExecuteWorkflow = () => {
		if (!json.isWorkflowValid) {
			return;
		}
		const result = runWorkflow(JsonDataState);
		JsonDataState.outputValue = result;
		JsonDataState.isExecutionDone = true;
	};

	const onResetWorkflow = () => {
		clearWorkflow();
	};

	const isExecuteWorkflowButtonDisabled = json.isWorkflowValid
		? 'bg-green-600 p-2 hover:bg-green-400'
		: 'bg-gray-400 p-2';

	return (
		<div className='flex flex-col my-1.5 gap-2 bg-slate-600 p-2 border-l-4 border-l-blue-500'>
			<p>
				<span className='font-bold italic'>Click</span> on the following actions
				below to execute the workflow:
			</p>
			<button
				className={isExecuteWorkflowButtonDisabled}
				onClick={onExecuteWorkflow}
				disabled={!json.isWorkflowValid}
			>
				Execute Workflow
			</button>
			<button
				className='bg-red-600 p-2 hover:bg-red-400'
				onClick={onResetWorkflow}
			>
				Reset Workflow
			</button>
		</div>
	);
};

export default WorkflowActions;
