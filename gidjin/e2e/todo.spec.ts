import { test, expect } from '@playwright/test';

test.beforeEach(async ({ page }) => {
  const url = process.env.TEST_URL || 'http://localhost:5173';
  await page.goto(url);
});

test('user can add todos', async ({ page }) => {
  await page.getByTestId('add-input').fill('buy milk');
  await page.getByTestId('add-button').click();
  await expect(page.getByTestId('item-0')).toHaveText('0. buy milkEditDelete');
  await page.getByTestId('add-input').fill('buy fruit');
  await page.getByTestId('add-button').click();
  await expect(page.getByTestId('item-1')).toHaveText('1. buy fruitEditDelete');
  await page.getByTestId('add-input').fill('buy bread');
  await page.getByTestId('add-button').click();
  await expect(page.getByTestId('item-2')).toHaveText('2. buy breadEditDelete');
});

test('user can edit todos', async ({ page }) => {
  await page.getByTestId('add-input').fill('buy milk');
  await page.getByTestId('add-button').click();
  await page.getByTestId('add-input').fill('buy fruit');
  await page.getByTestId('add-button').click();
  await expect(page.getByTestId('item-0')).toHaveText('0. buy milkEditDelete');
  await expect(page.getByTestId('edit-button-0')).toBeVisible();
  await page.getByTestId('edit-button-0').click();
  await expect(page.getByTestId('add-input')).toHaveValue('buy milk');
  await expect(page.getByTestId('add-button')).toHaveText('Save');
  await page.getByTestId('add-input').fill('buy oat milk');
  await page.getByTestId('add-button').click();
  await expect(page.getByTestId('item-0')).toHaveText('0. buy oat milkEditDelete');
});

test('user can delete todos', async ({ page }) => {
  await page.getByTestId('add-input').fill('buy milk');
  await page.getByTestId('add-button').click();
  await page.getByTestId('add-input').fill('buy fruit');
  await page.getByTestId('add-button').click();
  await expect(page.getByTestId('item-0')).toHaveText('0. buy milkEditDelete');
  await expect(page.getByTestId('item-1')).toHaveText('1. buy fruitEditDelete');
  await page.getByTestId('delete-button-0').click();
  await expect(page.getByTestId('item-0')).toHaveText('0. buy fruitEditDelete');
});
