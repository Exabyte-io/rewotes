import { Handle } from 'reactflow';

/***
 * This component is the Workflow End Node Component that is added onto the flowchart,
 * it implies that this is the end of the workflow and if the workflow is
 * executed and running, once the workflow reaches this component,
 * the workflow is expected to stop running, and end the execution safely.
 */
const WorkflowEndNode = (props) => {
	return (
		<>
			<Handle isConnectable position='top' type='target' />
			<div className='bg-green-600 border-slate-600 p-1 border-solid border-2 capitalize text-slate-700'>
				<p>{props.data.label}</p>
			</div>
		</>
	);
};

export default WorkflowEndNode;
