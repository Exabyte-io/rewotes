import React, { useState, useCallback, useRef, useEffect } from 'react';
import ReactFlow, {
  Controls,
  addEdge,
  type Node,
  type Edge,
  ReactFlowProvider,
  useNodesState,
  useEdgesState,
} from 'react-flow-renderer';
import DraggableItem from './DraggableItem';
import { BASE_CONDUCTIVITY, BATCH_SIZE, measurePConductivity, measureNConductivity, MARGIN } from './utils';
import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';
import { CanvasWrapper, DesignerContainer, JSONContainer, ResetButton, StyledAside } from './styles';

type DopantType = 'P-type' | 'N-type';

interface DesignerProps {
  dopingLevel: number;
  setDopingLevel: React.Dispatch<React.SetStateAction<number>>;
  desiredConductivity: number;
  setDesiredConductivity: React.Dispatch<React.SetStateAction<number>>;
  disabledPDoping: boolean;
  setDisablePDoping: React.Dispatch<React.SetStateAction<boolean>>;
  disabledNDoping: boolean;
  setDisableNdoping: React.Dispatch<React.SetStateAction<boolean>>;
  isFormSubmitted: boolean;
  setIsFormSubmitted: React.Dispatch<React.SetStateAction<boolean>>;
}

const initialNodes: Node[] = [];
const initialEdges: Edge[] = [];
const initialLoopData: any[] = [];

export const FlowchartDesigner = ({
  isFormSubmitted,
  setIsFormSubmitted,
  desiredConductivity,
  dopingLevel,
  setDopingLevel,
  disabledNDoping,
  disabledPDoping,
  setDisableNdoping,
  setDisablePDoping,
  setDesiredConductivity,
}: DesignerProps) => {
  const [nodes, setNodes, onNodesChange] = useNodesState(initialNodes);
  const [edges, setEdges, onEdgesChange] = useEdgesState(initialEdges);
  const [loopData, setLoopData] = useState<any[]>([]);
  const [feedback, setFeedback] = useState<any[]>([]);
  const reactFlowWrapper = useRef<HTMLDivElement>(null);
  const [batchSize, setBatchSize] = useState<number>(BATCH_SIZE);
  const [dopantType, setDopantType] = useState<DopantType | undefined>(undefined);
  const [goalAchieved, setGoalAchieved] = useState<boolean>(false);
  const [disableInc, setDisableInc] = useState<boolean>(false);
  const [disableDec, setDisableDec] = useState<boolean>(false);
  const [incAutoLoop, setIncAutoLoop] = useState<boolean>(false);
  const [decAutoLoop, setDecAutoLoop] = useState<boolean>(false);

  const initialConductivity = disabledNDoping
    ? measurePConductivity(dopingLevel)
    : disabledPDoping
    ? measureNConductivity(dopingLevel)
    : BASE_CONDUCTIVITY;
  const [conductivity, setConductivity] = useState<number>(initialConductivity);

  useEffect(() => {
    if (isFormSubmitted) {
      setFeedback((prev) => [...prev, { initialConductivity }, { desiredConductivity }]);
    }
  }, [isFormSubmitted, initialConductivity, desiredConductivity]);

  const isTypeDecisionAdded = (): boolean => {
    return nodes.some((node) => node.data.label === 'START_DOPING');
  };
  const isDopantTypeAdded = (): boolean => {
    return nodes.some((node) => node.data.label === 'P-TYPE_DOPING' || node.data.label === 'N-TYPE_DOPING');
  };
  const isCalcAdded = (): boolean => {
    return nodes.some((node) => node.data.label === 'CALCULATE-CONDUCTIVITY');
  };
  const isCheckAdded = (): boolean => {
    return nodes.some((node) => node.data.label === 'CHECK-CONDUCTIVITY');
  };
  const isIncrOrDecAdded = (): boolean => {
    return nodes.some((node) => node.data.label === 'INCREMENT' || node.data.label === 'DECREMENT');
  };
  const isBatchConditionAdded = (): boolean => {
    return nodes.some((node) => node.data.label === 'BATCH-CONDITION');
  };
  const isBatchChangeAdded = (): boolean => {
    return nodes.some(
      (node) =>
        node.data.label === 'DOUBLE_BATCH' || node.data.label === 'HALVE_BATCH' || node.data.label === 'MAINTAIN_BATCH',
    );
  };
  const isConclusion = (): boolean => {
    return nodes.some((node) => node.data.label === 'DONE');
  };

  const processNodeTypeLogic = (type: string) => {
    let additionalData = {};
    switch (type) {
      case 'START_DOPING':
        setFeedback((prev) => [...prev, 'Select doping type']);
        additionalData = {
          action: 'Select doping type',
        };
        break;
      case 'P-TYPE_DOPING': {
        const pType = 'P-type';
        setDopantType(pType);
        setFeedback((prev) => [...prev, 'Apply p type doping']);
        additionalData = {
          doping_type: pType,
          p_dopant: 'Boron',
          action: 'Apply p type doping',
        };
        break;
      }
      case 'N-TYPE_DOPING': {
        const nType = 'N-type';
        setDopantType(nType);
        setFeedback((prev) => [...prev, 'Apply n type doping']);
        additionalData = {
          doping_type: nType,
          n_dopant: 'Phosporus',
        };
        break;
      }
      case 'INCREMENT': {
        const highDopingLevel = dopingLevel + 10; // Simplified
        setDopingLevel(highDopingLevel);
        setFeedback((prev) => [...prev, `dopant increased from ${dopingLevel} to ${highDopingLevel}`]);
        additionalData = {
          dopantAdded: '10',
          oldDopingLevel: dopingLevel,
          newDopingLevel: highDopingLevel,
        };
        break;
      }
      case 'DECREMENT': {
        const lowDopingLevel = dopingLevel > 10 ? dopingLevel - 10 : 0; // Simplified
        setDopingLevel(lowDopingLevel);
        setFeedback((prev) => [...prev, `dopant decreased from ${dopingLevel} to ${lowDopingLevel}`]);
        additionalData = {
          dopantAdded: '0.1',
          oldDopingLevel: dopingLevel,
          newDopingLevel: lowDopingLevel,
        };
        break;
      }
      case 'CALCULATE-CONDUCTIVITY': {
        const measuredCond =
          dopantType === 'P-type' ? measurePConductivity(dopingLevel) : measureNConductivity(dopingLevel);
        setConductivity(measuredCond);
        setFeedback((prev) => [...prev, `The calculated conductivity is ${measuredCond}`]);
        if (Math.abs(measuredCond - desiredConductivity) <= MARGIN) {
          // Highlight the step in the flowchart, possibly changing its color or adding an icon.
          return;
        }
        additionalData = {
          action: 'Calculate conductiviy',
          measuredConductivity: measuredCond,
        };
        break;
      }
      case 'CHECK-CONDUCTIVITY': {
        const isDone = Math.abs(conductivity - desiredConductivity) <= MARGIN;
        if (isDone) {
          setGoalAchieved(true);
          setFeedback((prev) => [...prev, 'Desired conductivity range achieved!']);
          additionalData = {
            action: 'Check conductiviy',
            result: 'Desired conductivity range achieved!',
          };
        } else if (conductivity > desiredConductivity) {
          setDisableInc(true);
          setFeedback((prev) => [...prev, 'Conductivity too high. Decrease the dopant']);
          additionalData = {
            action: 'Check conductiviy',
            result: 'Conductivity too high. Decrease the dopant',
          };
        } else {
          setDisableDec(true);
          setFeedback((prev) => [...prev, 'Conductivity still low. Increase the dopant']);
          additionalData = {
            action: 'Check conductiviy',
            result: 'Conductivity still low. Increase the dopant',
          };
        }

        break;
      }
      case 'BATCH-CONDITION':
        setFeedback((prev) => [...prev, 'Select batch size adjustment']);
        additionalData = {
          action: 'Select batch size adjustment',
        };
        break;
      case 'DOUBLE_BATCH': {
        const doubledBatch = batchSize * 2;
        setBatchSize(doubledBatch);
        additionalData = {
          parameter: 'Increased batch size',
          oldBatchSize: BATCH_SIZE,
          newBatchSize: doubledBatch,
        };
        setFeedback((prev) => [...prev, `Batch size doubled from ${BATCH_SIZE} to ${doubledBatch}`]);
        break;
      }
      case 'HALVE_BATCH': {
        const halvedBatch = batchSize / 2;
        additionalData = {
          parameter: 'Decreased batch size',
          oldBatchSize: BATCH_SIZE,
          newValue: halvedBatch,
        };
        setFeedback((prev) => [...prev, `Batch size decreased from ${BATCH_SIZE} to ${halvedBatch}`]);
        break;
      }
      case 'MAINTAIN_BATCH':
        additionalData = {
          parameter: 'Maintained batch size',
          oldBatchSize: BATCH_SIZE,
          newValue: BATCH_SIZE,
        };
        setFeedback((prev) => [...prev, 'Maintained batch size']);
        break;
      case 'DONE':
        additionalData = {
          result: 'Ready for post-processing',
        };
        setFeedback((prev) => [...prev, 'Ready for post-processing']);
        break;
    }
    return additionalData;
  };

  const processNode = (type: string, event: React.DragEvent<HTMLDivElement>) => {
    const additionalData = processNodeTypeLogic(type);
    setNodes((prevNodes) => {
      if (type === 'START_DOPING' && prevNodes.some((node) => node.data.label === 'START_DOPING')) {
        toast.warning('You can only add one START_DOPING node.');
        return prevNodes;
      }

      if (
        (type === 'N-TYPE_DOPING' || type === 'P-TYPE_DOPING') &&
        (prevNodes.some((node) => node.data.label === 'N-TYPE_DOPING') ||
          prevNodes.some((node) => node.data.label === 'P-TYPE_DOPING'))
      ) {
        toast.warning('You can only add either (1) N-TYPE_DOPING or (1) P-TYPE_DOPING, not both.');
        return prevNodes;
      }

      const position = {
        x: event.clientX - (reactFlowWrapper.current?.offsetLeft || 0),
        y: event.clientY - (reactFlowWrapper.current?.offsetTop || 0),
      };

      const newNode: Node = {
        id: Math.random().toString(),
        type: 'custom',
        position,
        draggable: true,
        data: { label: type, ...additionalData },
      };

      if (prevNodes.length > 0) {
        const source = prevNodes[prevNodes.length - 1].id;
        const target = newNode.id;
        const newEdge: Edge = { id: `e-${source}-${target}`, source, target };
        setEdges((prevEdges) => [...prevEdges, newEdge]);

        // If the new node is "INCREMENT", create an edge to link it to the existing "CALC-CONDITION" node.
        if (type === 'INCREMENT') {
          const calcConditionNode = prevNodes.find((node) => node.data.label === 'CALCULATE-CONDUCTIVITY');
          if (calcConditionNode) {
            const newEdgeToCalc: Edge = {
              id: `e-${newNode.id}-${calcConditionNode.id}`,
              source: newNode.id,
              target: calcConditionNode.id,
            };
            setEdges((prevEdges) => [...prevEdges, newEdgeToCalc]);
            setIncAutoLoop(true);
          }
        }
        // If the new node is "DECREMENT", create an edge to link it to the existing "CALC-CONDITION" node.
        if (type === 'DECREMENT') {
          const calcConditionNode = prevNodes.find((node) => node.data.label === 'CALCULATE-CONDUCTIVITY');
          if (calcConditionNode) {
            const newEdgeToCalc: Edge = {
              id: `e-${newNode.id}-${calcConditionNode.id}`,
              source: newNode.id,
              target: calcConditionNode.id,
            };
            setEdges((prevEdges) => [...prevEdges, newEdgeToCalc]);
            setDecAutoLoop(true);
          }
        }
      }

      return [...prevNodes, newNode];
    });
  };

  const onDrop = useCallback(
    (event: React.DragEvent<HTMLDivElement>) => {
      event.preventDefault();

      const type = event.dataTransfer.getData('application/reactflow');
      processNode(type, event);
    },
    [reactFlowWrapper, nodes, conductivity],
  );

  const onDragOver = useCallback((event: React.DragEvent<HTMLDivElement>) => {
    event.preventDefault();
    event.dataTransfer.dropEffect = 'move';
  }, []);

  const executeIncLoop = () => {
    const calcLoopData = processNodeTypeLogic('CALCULATE-CONDUCTIVITY');
    setLoopData((prevData) => [...prevData, calcLoopData]);
    const checkLoopData = processNodeTypeLogic('CHECK-CONDUCTIVITY');
    setLoopData((prevData) => [...prevData, checkLoopData]);
    const incLoopData = processNodeTypeLogic('INCREMENT');
    setLoopData((prevData) => [...prevData, incLoopData]);
  };

  const executeDecLoop = () => {
    const calcLoopData = processNodeTypeLogic('CALCULATE-CONDUCTIVITY');
    setLoopData((prevData) => [...prevData, calcLoopData]);
    const checkLoopData = processNodeTypeLogic('CHECK-CONDUCTIVITY');
    setLoopData((prevData) => [...prevData, checkLoopData]);
    const incLoopData = processNodeTypeLogic('DECREMENT');
    setLoopData((prevData) => [...prevData, incLoopData]);
  };

  useEffect(() => {
    if (incAutoLoop) {
      const calcConditionNode = nodes.find((node) => node.data.label === 'CALCULATE-CONDUCTIVITY');
      if (calcConditionNode && Math.abs(conductivity - desiredConductivity) > MARGIN) {
        executeIncLoop();
      } else {
        setIncAutoLoop(false);
        setGoalAchieved(true);
        setFeedback((prev) => [...prev, 'Desired conducitivity range achieved!']);
      }
    }
  }, [incAutoLoop, nodes, conductivity]);

  useEffect(() => {
    if (decAutoLoop) {
      const calcConditionNode = nodes.find((node) => node.data.label === 'CALCULATE-CONDUCTIVITY');
      if (calcConditionNode && Math.abs(conductivity - desiredConductivity) > MARGIN) {
        executeDecLoop();
      } else {
        setDecAutoLoop(false);
        setGoalAchieved(true);
        setFeedback((prev) => [...prev, 'Desired conducitivity range achieved!']);
      }
    }
  }, [decAutoLoop, nodes, conductivity]);

  useEffect(() => {
    if ((goalAchieved && !incAutoLoop) || (goalAchieved && !decAutoLoop)) {
      // This is to prevent the alert
      // getting triggered too soon.
      setTimeout(() => {
        toast.success('Desired conductivity range achieved!', { position: toast.POSITION.TOP_CENTER });
      }, 600);
    }
  }, [goalAchieved, incAutoLoop, decAutoLoop]);

  useEffect(() => {
    if (disableDec) {
      toast.info('Conductivity is low. Increase the dopant', { position: toast.POSITION.TOP_CENTER });
    }
  }, [disableDec]);

  useEffect(() => {
    if (disableInc) {
      toast.info('Conductivity is too high. Decrease the dopant', { position: toast.POSITION.TOP_CENTER });
    }
  }, [disableInc]);

  const resetFlowchart = () => {
    setNodes(initialNodes);
    setEdges(initialEdges);
    setLoopData(initialLoopData);
    setFeedback(initialLoopData);
    setDopingLevel(0);
    setDopantType(undefined);
    setConductivity(0);
    setBatchSize(BATCH_SIZE);
    setGoalAchieved(false);
    setDisableInc(false);
    setDisableDec(false);
    setIncAutoLoop(false);
    setDecAutoLoop(false);
    setDisableNdoping(false);
    setDisablePDoping(false);
    setIsFormSubmitted(false);
    setDesiredConductivity(0);
  };

  return (
    <DesignerContainer>
      {/* Draggable items */}
      <h4 style={{ color: '#00435a', margin:'0' }}>Drag the nodes below onto the canvas to create your flowchart</h4>
      <StyledAside>
        <DraggableItem type="START_DOPING" isEnabled={isFormSubmitted && !isTypeDecisionAdded()} />
        <DraggableItem
          type="N-TYPE_DOPING"
          isEnabled={isTypeDecisionAdded() && !isDopantTypeAdded() && !disabledNDoping}
        />
        <DraggableItem
          type="P-TYPE_DOPING"
          isEnabled={isTypeDecisionAdded() && dopantType === undefined && !disabledPDoping}
        />
        <DraggableItem type="CALCULATE-CONDUCTIVITY" isEnabled={dopantType !== undefined && !isCalcAdded()} />
        <DraggableItem type="CHECK-CONDUCTIVITY" isEnabled={isCalcAdded() && !isCheckAdded()} />
        <DraggableItem
          type="INCREMENT"
          isEnabled={!goalAchieved && !disableInc && isCheckAdded() && !isIncrOrDecAdded()}
        />
        <DraggableItem
          type="DECREMENT"
          isEnabled={!goalAchieved && !disableDec && isCheckAdded() && !isIncrOrDecAdded()}
        />
        <DraggableItem type="BATCH-CONDITION" isEnabled={goalAchieved && !isBatchConditionAdded()} />
        <DraggableItem type="DOUBLE_BATCH" isEnabled={isBatchConditionAdded() && !isBatchChangeAdded()} />
        <DraggableItem type="HALVE_BATCH" isEnabled={isBatchConditionAdded() && !isBatchChangeAdded()} />
        <DraggableItem type="MAINTAIN_BATCH" isEnabled={isBatchConditionAdded() && !isBatchChangeAdded()} />
        <DraggableItem type="DONE" isEnabled={isBatchChangeAdded() && !isConclusion()} />
        <ResetButton onClick={resetFlowchart}>Reset Conditions & Flowchart</ResetButton>
      </StyledAside>
      <ToastContainer />
      <div style={{ display: 'flex' }}>
        <ReactFlowProvider>
          <CanvasWrapper ref={reactFlowWrapper}>
            <ReactFlow
              nodes={nodes}
              edges={edges}
              onDrop={onDrop}
              onDragOver={onDragOver}
              onNodesChange={onNodesChange}
              onEdgesChange={onEdgesChange}
              onConnect={(params) => {
                const newEdge = addEdge(params, edges);
                if (Array.isArray(newEdge)) {
                  setEdges((prevEdges) => [...prevEdges, ...newEdge]);
                } else {
                  setEdges((prevEdges) => [...prevEdges, newEdge]);
                }
              }}
            >
              <Controls />
            </ReactFlow>
          </CanvasWrapper>
        </ReactFlowProvider>
        <JSONContainer>
          <h3>Flowchart data (JSON)</h3>
          {JSON.stringify({ nodes, edges }, null, 2)}
        </JSONContainer>
        <JSONContainer>
          <h3>Feedback</h3>
          {JSON.stringify(feedback, null, 2)}
          <h3>Loop data</h3>
          {JSON.stringify(loopData, null, 2)}
        </JSONContainer>
      </div>
    </DesignerContainer>
  );
};

export default FlowchartDesigner;
