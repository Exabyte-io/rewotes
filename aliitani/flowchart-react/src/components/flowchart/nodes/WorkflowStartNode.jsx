import { Handle, Position } from 'reactflow';

const WorkflowStartNode = ({ data }) => {
	return (
		<>
			<div className='bg-green-600 border-slate-600 p-1 border-solid border-2 capitalize text-slate-700 text-center'>
				<p> {data.label} </p>
				<input
					type='number'
					placeholder='Enter initial value...'
					value={data.startValue}
					onChange={data.updateStartValue}
				/>
			</div>
			<Handle isConnectable position={Position.Bottom} type='source' />
		</>
	);
};

export default WorkflowStartNode;
