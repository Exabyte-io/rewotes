import React from 'react';
import { render, screen } from '@testing-library/react';
import JSONViewer from './JSONViewer';

describe('JSONViewer', () => {
  test('renders without crashing', () => {
    render(<JSONViewer nodes={[]} edges={[]} flows={[]} />);

    // Add more specific assertions based on your component's functionality
    // For example, you can check if certain elements are rendered, or if certain texts are displayed
  });
});