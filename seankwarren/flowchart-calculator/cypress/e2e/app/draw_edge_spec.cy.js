/// <reference types="cypress" />

describe('Draw Edge', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('draws the simplest flowchart with 1 edge', () => {
        // Click the clear button
        cy.get('.clear').click();
    
        // Place a new node on the flowchart
        cy.contains('button', 'in').dragAndDrop('.react-flow');
        cy.contains('button', 'out').dragAndDrop('.react-flow');
    
        // Connect the nodes
        cy.connectHandles(
            '.react-flow__nodes [data-testid=input-node] .handle.output',
            '.react-flow__nodes [data-testid=output-node] .handle.input'
        );
    
        // Check if the edge exists
        cy.get('.react-flow__edges')
            .find('.react-flow__edge')
            .should('have.length', 1);
    });
});