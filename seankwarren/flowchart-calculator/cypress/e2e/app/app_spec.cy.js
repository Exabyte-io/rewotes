/// <reference types="cypress" />

describe('App', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('user interaction elements renderer', () => {
        cy.get('button').contains('Clear');
        cy.get('button').contains('Save Flow');
        cy.get('button').contains('in');
        cy.get('button').contains('+');
        cy.get('button').contains('sin');
        cy.get('button').contains('>');
        cy.get('button').contains('out');
        cy.get('select').contains('Select a flow');
        cy.get('input');
    });
});