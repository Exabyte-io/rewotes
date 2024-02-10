import { describe, expect, test } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import App from "./App";

describe('App', () => {
  const user = userEvent.setup()

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
      await user.type(input, "test");
      await user.click(button);
      expect(screen.getByTestId("item-0")).toBeVisible();
      expect(screen.getByTestId("item-0")).toHaveTextContent("0. testEditDelete");
      await user.type(input, "second");
      await user.click(button);
      expect(screen.getByTestId("item-1")).toBeVisible();
      expect(screen.getByTestId("item-1")).toHaveTextContent("1. secondEditDelete");
    });
  });

  describe('deleting a todo', () => {
    test('renders delete button', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      await user.type(input, "test");
      await user.click(button);
      expect(screen.getByTestId("delete-button-0")).toBeVisible();
    });

    test('removes a todo', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      await user.type(input, "test");
      await user.click(button);
      expect(screen.getByTestId("delete-button-0")).toBeVisible();
      await user.click(screen.getByTestId("delete-button-0"));
      expect(screen.queryByTestId("item-0")).toBeNull();
    });
  });

  describe('editing a todo', () => {
    test('renders edit button', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      await user.type(input, "test");
      await user.click(button);
      expect(screen.getByTestId("edit-button-0")).toBeVisible();
    });

    test('edits a todo', async () => {
      render(<App />);
      const input = screen.getByTestId("add-input");
      const button = screen.getByTestId("add-button");
      await user.type(input, "test");
      await user.click(button);
      expect(screen.getByTestId("edit-button-0")).toBeVisible();
      await user.click(screen.getByTestId("edit-button-0"));

      expect(screen.getByTestId("add-input")).toHaveValue("test");
      expect(screen.getByTestId("add-button")).toHaveTextContent("Save");

      await user.clear(input);
      await user.type(input, "edited");
      await user.click(button);
      expect(screen.getByTestId("item-0")).toHaveTextContent("0. editedEditDelete");
      expect(screen.getByTestId("add-button")).toHaveTextContent("Add");
    });
  });
});

