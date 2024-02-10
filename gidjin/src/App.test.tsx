import { describe, expect, test } from 'vitest';
import { render, screen, fireEvent } from '@testing-library/react';
import App from "./App";

describe('App', () => {
  describe('adding a todo', () => {
    test('renders form', () => {
      render(<App />);
      expect(screen.getByTestId("add-input")).toBeVisible();
      expect(screen.getByTestId("add-button")).toBeVisible();
    });

    test('adds a todo', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      fireEvent.change(input, { target: { value: "test" } });
      await button.click();
      expect(screen.getByTestId("item-0")).toBeVisible();
      expect(screen.getByTestId("item-0").textContent).toEqual("0. testDelete");
    });
  });

  describe('deleting a todo', () => {
    test('renders delete button', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      fireEvent.change(input, { target: { value: "test" } });
      await button.click();
      expect(screen.getByTestId("delete-button-0")).toBeVisible();
    });

    test('removes a todo', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      fireEvent.change(input, { target: { value: "test" } });
      await button.click();
      expect(screen.getByTestId("delete-button-0")).toBeVisible();
      await screen.getByTestId("delete-button-0").click();
      expect(screen.queryByTestId("item-0")).toBeNull();
    });
  });
});

