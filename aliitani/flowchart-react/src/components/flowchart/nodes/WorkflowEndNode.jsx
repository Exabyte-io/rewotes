import { Handle } from 'reactflow';

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
