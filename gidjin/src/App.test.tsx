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
      expect(screen.getByTestId("item-0")).toHaveTextContent("0. testEditDelete");
      fireEvent.change(input, { target: { value: "second" } });
      await button.click();
      expect(screen.getByTestId("item-1")).toBeVisible();
      expect(screen.getByTestId("item-1")).toHaveTextContent("1. secondEditDelete");
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

  describe('editing a todo', () => {
    test('renders edit button', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      fireEvent.change(input, { target: { value: "test" } });
      await button.click();
      expect(screen.getByTestId("edit-button-0")).toBeVisible();
    });

    test('edits a todo', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      fireEvent.change(input, { target: { value: "test" } });
      await button.click();
      expect(screen.getByTestId("edit-button-0")).toBeVisible();
      await screen.getByTestId("edit-button-0").click();

      expect(screen.getByTestId("add-input")).toHaveValue("test");
      expect(screen.getByTestId("add-button")).toHaveTextContent("Save");

      fireEvent.change(input, { target: { value: "edited" } });
      await button.click();
      expect(screen.getByTestId("item-0")).toHaveTextContent("0. editedEditDelete");
      expect(screen.getByTestId("add-button")).toHaveTextContent("Add");
    });
  });
});

