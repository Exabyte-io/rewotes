import React from 'react';
import { render, fireEvent, screen } from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import SourceEditor from './SourceEditor';
import { SourceProvider } from '../../context/SourceContext';
import { AtomsProvider } from '../../context/AtomsContext';

const renderSourceEditor = () => {
  return render(
    <SourceProvider>
      <AtomsProvider>
        <SourceEditor />
      </AtomsProvider>
    </SourceProvider>
  );
};

describe('SourceEditor', () => {
  it('renders the component correctly', () => {
    renderSourceEditor();
    expect(screen.getByRole('heading', { level: 4 })).toHaveTextContent('Source editor');
  });

  it('displays defalut initial pattern validation message', () => {
    renderSourceEditor();
    expect(screen.getByText(/Wrong XYZ format/)).toBeInTheDocument();
  });

  it('updates the pattern validation message when the textarea content changes', () => {
    renderSourceEditor();
    const textArea = screen.getByLabelText('source-editor-textarea');

    fireEvent.change(textArea, {
      target: { value: '2\nComment line\nH 0 0 0\nH 0 0 1' },
    });
    expect(screen.getByText('Correct XYZ format')).toBeInTheDocument();

    // Comment line missing: wrong
    fireEvent.change(textArea, {
        target: { value: '2\nH 0 0 0\nH 0 0 1' }, 
      });
    expect(screen.getByText('Wrong XYZ format')).toBeInTheDocument();
    
    // NaN in atom lines: wrong
    fireEvent.change(textArea, {
        target: { value: '2\nComment line\nH mistake 0 0\nH 0 0 1' },
      });
    expect(screen.getByText('Wrong XYZ format')).toBeInTheDocument();
    
    // Atom lines number less then expected: wrong
    fireEvent.change(textArea, {
        target: { value: '3\nComment line\nH 0 0 0\nH 0 0 1' },
      });
    expect(screen.getByText('Wrong XYZ format')).toBeInTheDocument();
  });

  it('sets the correct Tag intent based on the validity of the XYZ format', () => {
    renderSourceEditor();

      // Default tag intent should be 'warning'

      const tag = screen.getByTestId('xyz-validity-tag');
    expect(tag).toHaveClass('bp4-intent-warning');
    expect(tag).toHaveTextContent('Wrong XYZ format');

    // Change the text area to a valid XYZ format
    fireEvent.change(screen.getByLabelText('source-editor-textarea'), {
      target: { value: '2\nComment line\nH 0 0 0\nH 0 0 1' },
    });

   // Tag intent should be 'success'
    expect(tag).toHaveClass('bp4-intent-success');
    expect(tag).toHaveTextContent('Correct XYZ format');
    });
});

