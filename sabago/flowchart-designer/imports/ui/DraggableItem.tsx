import React from 'react';
import { FlowchartElement as Element } from './styles';

const DraggableItem = ({ type, isEnabled }: { type: string; isEnabled: boolean }) => {
  let displayLabel = type;
  switch (type) {
    case 'START_DOPING':
      displayLabel = 'Doping Type?';
      break;
    case 'P-TYPE_DOPING':
      displayLabel = 'Apply (Boron) p-type Doping';
      break;
    case 'N-TYPE_DOPING':
      displayLabel = 'Apply (Phosphorus) n-type Doping';
      break;
    case 'INCREMENT':
      displayLabel = '+10 Doping';
      break;
    case 'DECREMENT':
      displayLabel = '-10 Doping';
      break;
    case 'DOUBLE_BATCH':
      displayLabel = 'Double Si Batch Size';
      break;
    case 'HALVE_BATCH':
      displayLabel = 'Halve Si Batch Size';
      break;
    case 'MAINTAIN_BATCH':
      displayLabel = 'Maintain Batch Size';
      break;
    case 'CALCULATE-CONDUCTIVITY':
      displayLabel = 'Calc Si Conductivity';
      break;
    case 'CHECK-CONDUCTIVITY':
      displayLabel = 'Check conductivity';
      break;
    case 'BATCH-CONDITION':
      displayLabel = 'Double or half the batch size?';
      break;
    case 'DONE':
      displayLabel = 'Ready for post processing';
      break;
    default:
      break;
  }

  return (
    <Element
      draggable={isEnabled}
      $isDisabled={isEnabled}
      onDragStart={(event) => {
        event.dataTransfer.setData('application/reactflow', type);
        event.dataTransfer.effectAllowed = 'move';
      }}
    >
      {displayLabel}
    </Element>
  );
};

export default DraggableItem;
