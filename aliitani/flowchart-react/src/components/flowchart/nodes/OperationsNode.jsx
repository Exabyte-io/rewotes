import { Handle } from 'reactflow';

const OperationsNode = ({ data }) => {
	return (
		<>
			<Handle isConnectable position='top' type='target' />
			<div className='bg-slate-500 border-slate-600 p-1 border-solid border-2 capitalize text-slate-400'>
				{data.label}{' '}
				<input
					type='number'
					placeholder='by this amount'
					onChange={data.updateOperationsValue}
				/>
			</div>
			<Handle isConnectable position='bottom' type='source' />
		</>
	);
};

export default OperationsNode;
