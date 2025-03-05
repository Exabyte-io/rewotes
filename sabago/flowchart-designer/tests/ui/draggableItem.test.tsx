import React from 'react';
import { render, screen, fireEvent } from '@testing-library/react';
import DraggableItem from '../../imports/ui/DraggableItem';

describe('DraggableItem Component', () => {
  it('renders without crashing', () => {
    render(<DraggableItem type="START_DOPING" isEnabled={true} />);
  });

  const testCases = [
    { type: 'START_DOPING', expectedLabel: 'Doping Type?' },
    { type: 'P-TYPE_DOPING', expectedLabel: 'Apply (Boron) p-type Doping' },
    { type: 'N-TYPE_DOPING', expectedLabel: 'Apply (Phosphorus) n-type Doping' },
    { type: 'INCREMENT', expectedLabel: '+10 Doping' },
    { type: 'DECREMENT', expectedLabel: '-10 Doping' },
    { type: 'DOUBLE_BATCH', expectedLabel: 'Double Si Batch Size' },
    { type: 'HALVE_BATCH', expectedLabel: 'Halve Si Batch Size' },
    { type: 'MAINTAIN_BATCH', expectedLabel: 'Maintain Batch Size' },
    { type: 'CALCULATE-CONDUCTIVITY', expectedLabel: 'Calc Si Conductivity' },
    { type: 'CHECK-CONDUCTIVITY', expectedLabel: 'Check conductivity' },
    { type: 'BATCH-CONDITION', expectedLabel: 'Double or half the batch size?' },
    { type: 'DONE', expectedLabel: 'Ready for post processing' },
  ];

  testCases.forEach(({ type, expectedLabel }) => {
    it(`displays the correct label for type: ${type}`, () => {
      render(<DraggableItem type={type} isEnabled={true} />);
      expect(screen.getByText(expectedLabel)).toBeInTheDocument();
    });
  });

  it('sets draggable attribute based on isEnabled prop', () => {
    const { rerender } = render(<DraggableItem type="START_DOPING" isEnabled={true} />);
    expect(screen.getByText('Doping Type?')).toHaveAttribute('draggable', 'true');

    rerender(<DraggableItem type="START_DOPING" isEnabled={false} />);
    expect(screen.getByText('Doping Type?')).toHaveAttribute('draggable', 'false');
  });

  it('sets the correct drag data on drag start', () => {
    render(<DraggableItem type="START_DOPING" isEnabled={true} />);
    const element = screen.getByText('Doping Type?');

    // Mock DataTransfer
    class MockDataTransfer {
      setData = jest.fn();
      effectAllowed: string | null = null;
    }

    const mockDataTransfer = new MockDataTransfer();

    // Fire the dragStart event
    fireEvent.dragStart(element, {
      dataTransfer: mockDataTransfer,
    });

    expect(mockDataTransfer.setData).toHaveBeenCalledWith('application/reactflow', 'START_DOPING');
    expect(mockDataTransfer.effectAllowed).toBe('move');
  });
});
