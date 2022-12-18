import { Handle, Position } from 'reactflow';
import { useState } from 'react';
import {
	addNewCondition,
	addNewOperation,
	getNode,
	updateObjectValueInJsonViewer,
} from '../../../state/AppState';

const WorkflowStartNode = ({ id, data }) => {
	const [startValue, setStartValue] = useState(data.startValue);
	const [isConnected, setIsConnected] = useState(false);

	const onChangeHandler = (e) => {
		const valueToUpdate = Number(e.target.value);
		setStartValue(valueToUpdate);

		let targetNode = getNode(id);
		targetNode.data.startValue = valueToUpdate;
		updateObjectValueInJsonViewer(valueToUpdate, 'startValue');
	};

	const onConnectHandler = (params) => {
		const targetNode = getNode(params.target);
		if (isConnected) {
			return;
		}

		if (targetNode.type === 'operations') {
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
				targetNode.data.value,
				targetNode.type
			);
		}
		setIsConnected(true);
	};

	return (
		<>
			<div className='bg-green-600 border-slate-600 p-1 border-solid border-2 capitalize text-slate-700 text-center'>
				<p> {data.label} </p>
				<input
					className='text-right'
					type='number'
					placeholder='Enter initial value...'
					value={startValue}
					onChange={onChangeHandler}
				/>
			</div>
			<Handle
				isConnectable
				position={Position.Bottom}
				type='source'
				onConnect={onConnectHandler}
			/>
		</>
	);
};

export default WorkflowStartNode;
