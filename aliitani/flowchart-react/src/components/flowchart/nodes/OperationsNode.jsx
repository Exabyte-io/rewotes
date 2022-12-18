import { useState } from 'react';
import { Handle } from 'reactflow';
import {
	getNode,
	updateOperationByNode,
	addNewOperation,
	addNewCondition,
} from '../../../state/AppState';

const OperationsNode = ({ id, data }) => {
	const [operation, setOperation] = useState({ operation: '', value: 0 });
	const [isConnected, setIsConnected] = useState(false);

	const onChangeHandler = (e) => {
		const value = Number(e.target.value);
		setOperation({ operation: data.label, value: value });
		const targetNode = getNode(id);
		targetNode.data.value = value;
		updateOperationByNode(targetNode);
	};

	const onConnectHandler = (params) => {
		const targetNode = getNode(params.target);
		if (isConnected) {
			return;
		}

		console.log(targetNode.type);
		if (targetNode.type === 'operations') {
			// add into operations
			addNewOperation(
				id,
				targetNode.id,
				targetNode.data.operationType,
				targetNode.data.value
			);
		} else if (targetNode.type === 'conditionals') {
			addNewCondition(
				id,
				targetNode.id,
				targetNode.data.condition,
				targetNode.data.value
			);
		} else if (targetNode.type === 'endNode') {
			addNewOperation(id, targetNode.id, 'end_workflow', targetNode.data.value);
		}
		setIsConnected(true);
	};

	return (
		<>
			<Handle isConnectable position='top' type='target' />
			<div className='bg-slate-500 border-slate-600 p-1 border-solid border-2 capitalize text-slate-400'>
				{data.operationType}{' '}
				<input
					type='number'
					placeholder='by this amount'
					value={operation.value}
					onChange={onChangeHandler}
				/>
			</div>
			<Handle
				isConnectable
				position='bottom'
				type='source'
				onConnect={(params) => onConnectHandler(params)}
			/>
		</>
	);
};

export default OperationsNode;
