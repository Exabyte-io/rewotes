import { Handle, Position } from 'reactflow';

const ConditionalsNode = ({ data }) => {
	return (
		<div className='border-2 border-slate-600'>
			<div className='text-slate-700 bg-orange-400  p-2'>
				<Handle position={Position.Top} type='target' />
				<p>If the value is</p>
				<div className='flex flex-row flex-nowrap gap-2'>
					<p>is</p>
					<select
						className='bg-slate-500 text-white'
						onChange={data.updateConditionalValue}
					>
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
						onChange={data.updateExpectedValue}
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
				id='true'
				style={{ left: 10, background: '#555' }}
			/>
			<Handle
				type='source'
				position={Position.Bottom}
				id='false'
				style={{ right: 10, left: 'auto', background: '#555' }}
			/>
		</div>
	);
};

export default ConditionalsNode;
