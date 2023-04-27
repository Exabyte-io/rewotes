import { fireEvent, render, screen } from '@testing-library/react';
import StructureViewer from './StructureViewer';
import { AtomsProvider } from '../../context/AtomsContext';
import { SettingsProvider } from '../../context/SettingsContext';
import { SourceProvider } from '../../context/SourceContext';
import React from 'react';

jest.mock('@react-three/fiber');

describe('StructureViewer', () => {
  const renderStructureViewer = () => {
    return render(
      <AtomsProvider>
        <SettingsProvider>
          <SourceProvider>
            <StructureViewer />
          </SourceProvider>
        </SettingsProvider>
      </AtomsProvider>
    );
  };

  it('calls onAtomMove when dragging an atom', () => {
    const setSource = jest.fn();
    jest.spyOn(React, 'useContext').mockImplementation(() => ({
      source: '2\nComment\nH 0 0 0\nH 0 0 1',
      setSource,
    }));

    const draggableSphere = screen.getByTitle('atom-0');
    fireEvent.mouseDown(draggableSphere);
    fireEvent.mouseMove(draggableSphere, { clientX: 100, clientY: 100 });
    fireEvent.mouseUp(draggableSphere);

    expect(setSource).toHaveBeenCalled();
  });
});