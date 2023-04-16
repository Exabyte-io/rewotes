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
        cy.contains('button', '+').dragAndDrop('.react-flow');
        cy.contains('button', 'out').dragAndDrop('.react-flow');
    
        // Connect the nodes
        cy.connectHandles(
            '.react-flow__nodes [data-testid=input-node] .handle.source',
            '.react-flow__nodes [data-testid=binary-node] .handle.target.top'
        );

        cy.connectHandles(
            '.react-flow__nodes [data-testid=input-node] .handle.source',
            '.react-flow__nodes [data-testid=binary-node] .handle.target.bottom'
        );

        cy.connectHandles(
            '.react-flow__nodes [data-testid=binary-node] .handle.source',
            '.react-flow__nodes [data-testid=output-node] .handle.target'
        );
    
        cy.get('[data-testid=input-node] input')
            .clear({ force: true }) // Clear the current value of the input element
            .type('42', { force: true }); // Set the new value to '42'
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '84')
    });
});