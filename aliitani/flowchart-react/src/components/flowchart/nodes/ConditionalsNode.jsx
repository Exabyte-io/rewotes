import { useState } from 'react';
import { Handle, Position } from 'reactflow';
import {
	getNode,
	updateConditionByNode,
	updateTruthStateConditionByNode,
	addNewOperation,
} from '../../../state/AppState';

/***
 * This component is the Conditional Component that is added onto the flowchart, it contains the conditional operators
 * that run conditional checks on the workflow value. The conditional operators are the basic conditional checks.
 * It contains two output connections, one for the truthState and the other is for the falseState.
 */
const ConditionalsNode = ({ id, data }) => {
	const [isConnected, setIsConnected] = useState({
		truthState: false,
		falseState: false,
	});

	const onChangeSelectOptionHandler = (e) => {
		const targetNode = getNode(id);
		const valueToChange = e.target.value;

		targetNode.data.condition = valueToChange;
		updateConditionByNode(targetNode, 'condition');
	};

	const onChangeValueHandler = (e) => {
		const targetNode = getNode(id);
		const valueToChange = Number(e.target.value);

		targetNode.data.value = valueToChange;
		updateConditionByNode(targetNode, 'value');
	};

	const onConnectHandler = (params) => {
		if (isConnected[params.sourceHandle]) {
			return;
		}

		const sourceNode = getNode(params.source);
		const targetNode = getNode(params.target);

		if (
			sourceNode.type === 'conditionals' &&
			targetNode.type === 'conditionals'
		) {
			return;
		}

		console.log(params, targetNode);

		if (targetNode.type === 'operations') {
			addNewOperation(
				params.source,
				targetNode.id,
				targetNode.data.operationType,
				targetNode.data.value
			);
		}

		updateTruthStateConditionByNode(
			params.source,
			params.sourceHandle,
			targetNode.id
		);

		setIsConnected((c) => ({ ...c, [params.sourceHandle]: true }));
	};

	return (
		<div className='border-2 border-slate-600'>
			<div className='text-slate-700 bg-orange-400  p-2'>
				<Handle position={Position.Top} type='target' />
				<p>If the value is</p>
				<div className='flex flex-row flex-nowrap gap-2'>
					<p>is</p>
					<select
						className='bg-slate-500 text-white'
						onChange={onChangeSelectOptionHandler}
						defaultValue='<'
					>
						<option value='<' disabled>
							Select Value
						</option>
						<option value='<'>{'< - Less than'}</option>
						<option value='>'>{'> - Greater than'}</option>
						<option value='<='>{'<= - Less than/Equal to'}</option>
						<option value='>='>{'>= - Greater than/Equal to'}</option>
						<option value='=='>{'== - Equal to'}</option>
						<option value='!='>{'!= - Not Equal to'}</option>
					</select>
				</div>
				<div className='flex flex-col flex-wrap'>
					<label htmlFor='expectedValue'>the expected value, which is</label>
					<input
						type='number'
						id='expectedValue'
						defaultValue='0'
						onChange={onChangeValueHandler}
					/>
				</div>
			</div>
			<div className='flex flex-row flex-nowrap w-full bg-slate-400 px-2 justify-between'>
				<p>True</p>
				<p>False</p>
			</div>
			<Handle
				type='source'
				position={Position.Bottom}
				id='truthState'
				style={{ left: 10, background: '#555' }}
				onConnect={onConnectHandler}
			/>
			<Handle
				type='source'
				position={Position.Bottom}
				id='falseState'
				style={{ right: 10, left: 'auto', background: '#555' }}
				onConnect={onConnectHandler}
			/>
		</div>
	);
};

export default ConditionalsNode;
