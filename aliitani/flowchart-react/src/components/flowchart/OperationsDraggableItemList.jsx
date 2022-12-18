import React from 'react';

const OPERATIONS = {
	INCREMENT: 'increment',
	DECREMENT: 'decrement',
	MULTIPLY: 'multiply',
	DIVIDE: 'divide',
	CONDITIONALS: 'conditionals',
};
const DragOperation = ({ label }) => {
	const onDragStart = (event) => {
		const type = label === 'conditionals' ? 'conditionals' : 'operations';
		event.dataTransfer.setData('application/reactflow', type);
		event.dataTransfer.setData('label', label);
		event.dataTransfer.effectAllowed = 'move';
	};

	const classNameOperation =
		label === OPERATIONS.CONDITIONALS
			? 'bg-orange-500 p-2 hover:bg-orange-400 cursor-grab capitalize'
			: 'bg-slate-500 p-2 hover:bg-slate-400 cursor-grab capitalize';
	return (
		<div
			className={classNameOperation}
			draggable
			onDragStart={(event) => onDragStart(event)}
		>
			{label}
		</div>
	);
};

const OperationsDraggableItemList = () => {
	return (
		<div className='flex flex-col flex-nowrap gap-2 bg-slate-600 p-2 border-l-4 border-l-blue-500'>
			<p>
				<span className='font-bold italic'>Drag</span> any of these operations
				on the screen to include in the workflow:
			</p>
			<DragOperation label={OPERATIONS.INCREMENT} />
			<DragOperation label={OPERATIONS.DECREMENT} />
			<DragOperation label={OPERATIONS.MULTIPLY} />
			<DragOperation label={OPERATIONS.DIVIDE} />
			<DragOperation label={OPERATIONS.CONDITIONALS} />
		</div>
	);
};

export default OperationsDraggableItemList;
