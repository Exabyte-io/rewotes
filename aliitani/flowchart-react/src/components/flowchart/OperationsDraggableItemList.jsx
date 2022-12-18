import React from 'react';

const OPERATIONS = {
	INCREMENT: 'increment',
	DECREMENT: 'decrement',
	MULTIPLY: 'multiply',
	DIVIDE: 'divide',
	CONDITIONALS: 'conditionals',
};
const DragOperation = ({ nodeType }) => {
	const onDragStart = (event, nodeType) => {
		event.dataTransfer.setData('application/reactflow', nodeType);
		event.dataTransfer.effectAllowed = 'move';
	};

	const classNameOperation =
		nodeType === OPERATIONS.CONDITIONALS
			? 'bg-orange-500 p-2 hover:bg-orange-400 cursor-grab capitalize'
			: 'bg-slate-500 p-2 hover:bg-slate-400 cursor-grab capitalize';
	return (
		<div
			className={classNameOperation}
			draggable
			onDragStart={(event) => onDragStart(event, nodeType)}
		>
			{nodeType}
		</div>
	);
};

const OperationsDraggableItemList = () => {
	return (
		<div className='flex flex-col flex-nowrap gap-2'>
			<div className='flex flex-col flex-nowrap gap-2 bg-slate-600 p-2 border-l-4 border-l-blue-500'>
				<p>
					<span className='font-bold italic'>Drag</span> any of these operations
					on the screen to include in the workflow:
				</p>
				<DragOperation nodeType={OPERATIONS.INCREMENT} />
				<DragOperation nodeType={OPERATIONS.DECREMENT} />
				<DragOperation nodeType={OPERATIONS.MULTIPLY} />
				<DragOperation nodeType={OPERATIONS.DIVIDE} />
				<DragOperation nodeType={OPERATIONS.CONDITIONALS} />
			</div>
		</div>
	);
};

export default OperationsDraggableItemList;
