/// <reference types="cypress" />

describe('Save Flow', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('saves a flow to MongoDB and then reloads it', () => {
        // Click the clear button
        cy.get('.clear').click();
        
        // Place a new node on the flowchart
        cy.contains('button', 'in').dragAndDrop('.react-flow');

        // Type a flow name into the input field
        cy.get('.flowname').type("automated save test")

        // Click the save button
        cy.get('.saveflow').click();

        // Click the clear button
        cy.get('.clear').click();

        // Deselect a flow from the dropdown selection
        cy.get('.flow-dropdown').select('Select a flow')

        // Select the new flow from the dropdown selection
        cy.get('.flow-dropdown').select('automated save test')

        cy.get('.react-flow__nodes')
            .find('[data-testid=input-node]')
            .should('exist');
    });
});