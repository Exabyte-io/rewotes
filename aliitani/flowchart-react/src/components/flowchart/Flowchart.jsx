import React, { useCallback, useState, useRef, useEffect } from 'react';
import ReactFlow, {
	Controls,
	Background,
	ReactFlowProvider,
	useNodesState,
	useEdgesState,
	MiniMap,
	addEdge,
} from 'reactflow';
import 'reactflow/dist/style.css';
import './flowchart.css';

import Notes from '../Notes';
import WorkflowActions from './WorkflowActions';
import ConditionalsNode from './nodes/ConditionalsNode';
import WorkflowStartNode from './nodes/WorkflowStartNode';
import WorkflowEndNode from './nodes/WorkflowEndNode';
import OperationsNode from './nodes/OperationsNode';
import OperationsDraggableItemList from './OperationsDraggableItemList';
import { subscribe, useSnapshot } from 'valtio';
import {
	Nodes,
	addNewNode,
	getNode,
	updateNodePosition,
	updateObjectValueInJsonViewer,
} from '../../state/AppState';

const reactFlowStyle = {
	background: '#cecece',
	height: '100%',
	width: '100%',
};

const nodeTypes = {
	operations: OperationsNode,
	conditionals: ConditionalsNode,
	startNode: WorkflowStartNode,
	endNode: WorkflowEndNode,
};

/***
 * This component holds the entire flowchart ui. It handles events triggered against the flowchart.
 */
const Flowchart = () => {
	/** The connections state is added to only allow one connection going out of a node from source -> target. */
	const [connections, setConnections] = useState({});

	const nodeState = useSnapshot(Nodes);

	const [nodes, setNodes] = useNodesState(Object.values(nodeState));
	const [edges, setEdges] = useEdgesState([]);

	const reactFlowWrapper = useRef(null);
	const [reactFlowInstance, setReactFlowInstance] = useState(null);

	const unsubscribeNodeState = subscribe(Nodes, () => {
		setNodes(Object.values(Nodes));
	});
	useEffect(() => () => unsubscribeNodeState(), []);

	const onDragOverHandler = useCallback((event) => {
		event.preventDefault();
		event.dataTransfer.dropEffect = 'move';
	}, []);

	const onConnectHandler = (params) => {
		const node = getNode(params.source);
		const isConditional = node.type === 'conditionals';

		if (
			isConditional &&
			connections[params.source] &&
			connections[params.source][params.sourceHandle]
		) {
			console.log('inside');
			return;
		}

		if (!isConditional && connections[params.source]) {
			return;
		}

		if (params.source === 'node_start_value') {
			updateObjectValueInJsonViewer(true, 'startNodeConnected');
		}

		if (params.target === 'node_result_value') {
			updateObjectValueInJsonViewer(true, 'workflowValid');
		}

		const sourceNode = getNode(params.source);
		const targetNode = getNode(params.target);

		if (
			sourceNode.type === 'conditionals' &&
			targetNode.type === 'conditionals'
		) {
			return;
		}

		if (isConditional) {
			setConnections((c) =>
				Object.assign(c, {
					[params.source]: { ...c[params.source], [params.sourceHandle]: true },
				})
			);
		} else {
			setConnections((c) => Object.assign(c, { [params.source]: true }));
		}
		setEdges((eds) => addEdge(params, eds));
	};

	const onDropHandler = useCallback(
		(event) => {
			event.preventDefault();
			const reactFlowBounds = reactFlowWrapper.current.getBoundingClientRect();
			const type = event.dataTransfer.getData('application/reactflow');
			const label = event.dataTransfer.getData('label');

			if (typeof type === 'undefined' || !type) {
				return;
			}

			addNewNode(
				type,
				reactFlowInstance.project({
					x: event.clientX - reactFlowBounds.left,
					y: event.clientY - reactFlowBounds.top,
				}),
				{
					label: `${type}`,
					operationType: label,
					value: 0,
				}
			);
		},
		[reactFlowInstance]
	);

	const onNodesChangeHandler = useCallback((changes) => {
		const [nodeToUpdate] = changes;
		if (nodeToUpdate.type !== 'position') {
			return;
		}
		updateNodePosition(nodeToUpdate.position, nodeToUpdate.id);
	}, []);

	return (
		<div>
			<ReactFlowProvider>
				<div id='flowchart-container' ref={reactFlowWrapper}>
					<ReactFlow
						style={reactFlowStyle}
						fitView
						onInit={setReactFlowInstance}
						onConnect={onConnectHandler}
						onDrop={onDropHandler}
						onDragOver={onDragOverHandler}
						onNodesChange={onNodesChangeHandler}
						nodes={nodes}
						nodeTypes={nodeTypes}
						edges={edges}
					>
						<Background variant='lines' />
						<Controls />
						<MiniMap />
					</ReactFlow>
				</div>
				<div className='flex flex-col flex-nowrap gap-2'>
					<Notes />
					<OperationsDraggableItemList />
					<WorkflowActions />
				</div>
			</ReactFlowProvider>
		</div>
	);
};

export default Flowchart;
