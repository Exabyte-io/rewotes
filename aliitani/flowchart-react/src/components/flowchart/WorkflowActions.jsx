import { clearWorkflow } from '../../state/AppState';

/***
 * This component holds the workflow actions that gets triggered against the flowchart.
 * - ExecuteWorkflow - executed the workflow based on the nodes(elements) rendered into the flowchart.
 * - ResetWorkflow - resets the workflow and clears out any state.
 */
const WorkflowActions = () => {
	const onExecuteWorkflow = () => {
		console.log('executing workflow...');
	};

	const onResetWorkflow = () => {
		clearWorkflow();
	};

	return (
		<div className='flex flex-col flex-grow gap-2 bg-slate-600 p-2 border-l-4 border-l-blue-500'>
			<p>
				<span className='font-bold italic'>Click</span> on the following actions
				below to execute the workflow:
			</p>
			<button
				className='bg-green-600 p-2 hover:bg-green-400'
				onClick={onExecuteWorkflow}
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
