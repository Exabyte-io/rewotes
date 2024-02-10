import { describe, expect, test } from 'vitest';
import { render, screen } from '@testing-library/react';
import App from "./App";


describe('App', () => {
  test('renders headline', () => {
    render(<App />);
    const headline = screen.getByText(/Vite \+ React/i);
    expect(headline).toBeInTheDocument();
  });
});

