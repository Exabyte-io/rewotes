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

import WorkflowActions from './WorkflowActions';
import ConditionalsNode from './nodes/ConditionalsNode';
import WorkflowStartNode from './nodes/WorkflowStartNode';
import WorkflowEndNode from './nodes/WorkflowEndNode';
import OperationsNode from './nodes/OperationsNode';
import OperationsDraggableItemList from './OperationsDraggableItemList';

const reactFlowStyle = {
	background: '#cecece',
	height: '100%',
	width: '100%',
};

const nodeTypes = {
	increment: OperationsNode,
	decrement: OperationsNode,
	multiply: OperationsNode,
	divide: OperationsNode,
	conditionals: ConditionalsNode,
	startNode: WorkflowStartNode,
	endNode: WorkflowEndNode,
};

let id = 0;
const getId = () => `node_${id++}`;

const Flowchart = ({ updateJsonData, jsonData, resetJson }) => {
	const [operationalNodeValueChanges, setOperationalNodeValueChanges] =
		useState({});
	const [conditionalNodeValueChanges, setConditionalNodeValueChanges] =
		useState({});
	const [startValue, setStartValue] = useState(0);
	const [outputValue, setOutputValue] = useState(0);

	const reactFlowWrapper = useRef(null);
	const [edges, setEdges, _] = useEdgesState([]);
	const [reactFlowInstance, setReactFlowInstance] = useState(null);
	const [nodes, setNodes, onNodesChange] = useNodesState([
		{
			id: 'node_start_value',
			type: 'startNode',
			data: {
				label: 'Start',
				updateStartValue: (e) => updateStartValueEventHandler(e),
				updateExpectedValue: (e) =>
					updateExpectedValueEventHandler(e, newNodeId),
				updateConditionalValue: (e) =>
					updateConditionalValueEventHandler(e, newNodeId),
				updateOperationsValue: (e) =>
					updateOperationsValueEventHandler(e, newNodeId),
			},
			position: { x: 250, y: 5 },
		},
		{
			id: 'node_result_value',
			type: 'endNode',
			data: {
				label: 'End of Workflow',
			},
			position: { x: 250, y: 150 },
		},
	]);

	const clearNodes = () => {
		setNodes([
			{
				id: 'node_start_value',
				type: 'startNode',
				data: {
					label: 'Start',
					startValue: 0,
					updateStartValue: (e) => updateStartValueEventHandler(e),
					updateExpectedValue: (e) =>
						updateExpectedValueEventHandler(e, newNodeId),
					updateConditionalValue: (e) =>
						updateConditionalValueEventHandler(e, newNodeId),
					updateOperationsValue: (e) =>
						updateOperationsValueEventHandler(e, newNodeId),
				},
				position: { x: 250, y: 5 },
			},
			{
				id: 'node_result_value',
				type: 'endNode',
				data: {
					label: 'End of Workflow',
				},
				position: { x: 250, y: 150 },
			},
		]);
		setEdges([]);
		resetJson();
		setStartValue(0);
		setOutputValue(0);
		setConditionalNodeValueChanges({});
		setOperationalNodeValueChanges({});
	};

	const onDragOverHandler = useCallback((event) => {
		event.preventDefault();
		event.dataTransfer.dropEffect = 'move';
	}, []);

	const onConnectHandler = useCallback(
		(params) => {
			setEdges((eds) => addEdge(params, eds));
		},
		[reactFlowInstance, nodes]
	);

	const updateStartValueEventHandler = useCallback(
		(e) => {
			e.preventDefault();
			setStartValue(e.target.value);

			setNodes((node) => {
				return node.map((n) => {
					if (n.id === 'node_start_value') {
						return {
							...n,
							data: {
								...n.data,
								startValue: e.target.value,
							},
						};
					}

					return { ...n };
				});
			});
		},
		[reactFlowInstance, nodes]
	);

	const updateOperationsValueEventHandler = useCallback(
		(e, nodeId) => {
			e.preventDefault();
			const operations = jsonData.operations;
			const hasTargetNode = operations.filter(
				(operation) => operation.target === nodeId
			);
			console.log(operations, hasTargetNode);
			setOperationalNodeValueChanges((node) => ({
				...node,
				[nodeId]: {
					...operationalNodeValueChanges[nodeId],
					value: e.target.value,
					nodeId: nodeId,
				},
			}));
		},
		[reactFlowInstance, nodes, setNodes]
	);

	const updateExpectedValueEventHandler = useCallback(
		(e, nodeId) => {
			e.preventDefault();

			setConditionalNodeValueChanges((node) => ({
				...node,
				[nodeId]: {
					...conditionalNodeValueChanges[nodeId],
					expectedValue: e.target.value,
				},
			}));
		},
		[reactFlowInstance, nodes]
	);

	const updateConditionalValueEventHandler = useCallback(
		(e, nodeId) => {
			e.preventDefault();

			setConditionalNodeValueChanges((node) => ({
				...node,
				[nodeId]: {
					...conditionalNodeValueChanges[nodeId],
					conditionalValue: e.target.value,
				},
			}));
		},
		[reactFlowInstance, nodes]
	);

	const onDropHandler = useCallback(
		(event) => {
			event.preventDefault();

			const reactFlowBounds = reactFlowWrapper.current.getBoundingClientRect();
			const type = event.dataTransfer.getData('application/reactflow');
			if (typeof type === 'undefined' || !type) {
				return;
			}

			const position = reactFlowInstance.project({
				x: event.clientX - reactFlowBounds.left,
				y: event.clientY - reactFlowBounds.top,
			});

			const newNodeId = getId();
			const newNode = {
				id: newNodeId,
				type,
				position,
				data: {
					label: `${type}`,
					updateStartValue: (e) => updateStartValueEventHandler(e),
					updateExpectedValue: (e) =>
						updateExpectedValueEventHandler(e, newNodeId),
					updateConditionalValue: (e) =>
						updateConditionalValueEventHandler(e, newNodeId),
					updateOperationsValue: (e) =>
						updateOperationsValueEventHandler(e, newNodeId),
				},
			};

			setNodes((nds) => [...nds, newNode]);
		},
		[reactFlowInstance, nodes, setNodes]
	);

	useEffect(() => {
		updateJsonData((json) => ({ ...json, startValue: Number(startValue) }));
	}, [startValue]);

	useEffect(() => {
		if (!edges.length) {
			return;
		}

		const newEdge = edges[edges.length - 1];
		const nodeTarget = nodes.find((node) => node.id === newEdge.target);

		updateJsonData((json) => ({
			...json,
			operations: [
				...json.operations,
				{
					source: newEdge.source,
					target: newEdge.target,
					operation: nodeTarget.type,
					value: '',
				},
			],
		}));

		updateOperationByEdges();
	}, [edges]);

	useEffect(() => {
		console.log(
			nodes,
			conditionalNodeValueChanges,
			operationalNodeValueChanges
		);
	}, [nodes, conditionalNodeValueChanges, operationalNodeValueChanges]);

	useEffect(() => {
		updateOperationByEdges();
	}, [operationalNodeValueChanges]);

	const updateOperationByEdges = useCallback(() => {
		updateJsonData((json) => {
			const newOperations = json.operations.map((operation) => ({
				...operation,
				value: operationalNodeValueChanges[operation.target]?.value,
			}));

			console.log(newOperations);
			return {
				...json,
				operations: newOperations,
			};
		});
	}, [operationalNodeValueChanges]);

	const execute = () => {
		console.log('execute data', nodes, edges);
	};

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
						nodes={nodes}
						nodeTypes={nodeTypes}
						onNodesChange={onNodesChange}
						edges={edges}
					>
						<Background variant='lines' />
						<Controls />
						<MiniMap />
					</ReactFlow>
				</div>
				<div className='flex flex-col flex-nowrap gap-2'>
					<OperationsDraggableItemList />
					<WorkflowActions
						deleteWorkflow={clearNodes}
						executeWorkflow={execute}
					/>
				</div>
			</ReactFlowProvider>
		</div>
	);
};

export default Flowchart;
