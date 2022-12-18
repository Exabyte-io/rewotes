import { useState } from 'react';
import { Handle, Position } from 'reactflow';

const StartNode = () => {
	const [startValue, setStartValue] = useState(0);

	const onChangeHandler = (e) => {
		e.preventDefault();
		setStartValue(e.target.value);
	};

	return (
		<>
			<input type='text' onChange={onChangeHandler} />
			<Handle
				type='source'
				position={Position.Bottom}
				className='bg-green-500 text-slate-600 border-2 border-slate-400 p-2'
				onConnect={(params) => transmit(params)}
			/>
		</>
	);
};
