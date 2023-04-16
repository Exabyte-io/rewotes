/// <reference types="cypress" />

describe('Dark Mode Switch', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('toggles dark mode on and off', () => {
        // Check if the 'dark-mode' class is not applied initially
        cy.get('.react-flow').should('have.class', 'dark-mode');
        cy.get('.buttons-panel').should('have.class', 'dark-mode');
        cy.get('.json-viewer').should('have.class', 'dark-mode');
    
        // Toggle the dark mode switch
        cy.get('.darkmode-switch').click();
    
        // Check if the 'dark-mode' class is applied to the body
        cy.get('.react-flow').should('not.have.class', 'dark-mode');
        cy.get('.buttons-panel').should('not.have.class', 'dark-mode');
        cy.get('.json-viewer').should('not.have.class', 'dark-mode');
    
        // Toggle the dark mode switch back to off
        cy.get('.darkmode-switch').click();
    
        // Check if the 'dark-mode' class is removed from the body
        cy.get('.react-flow').should('have.class', 'dark-mode');
        cy.get('.buttons-panel').should('have.class', 'dark-mode');
        cy.get('.json-viewer').should('have.class', 'dark-mode');
    });
});