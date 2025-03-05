import React from 'react';
import { render, screen } from '@testing-library/react';
import NodeButtons from './NodeButtons';

describe('NodeButtons', () => {
  test('renders without crashing', () => {
    render(<NodeButtons handleDragStart={() => {}} isDarkMode={true} clearFlowchart={() => {}} clearFlows={() => {}} />);

    // Add more specific assertions based on your component's functionality
    // For example, you can check if certain elements are rendered, or if certain texts are displayed
  });
});