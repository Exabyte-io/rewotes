/// <reference types="cypress" />

describe('Basic user interaction', () => {
    beforeEach(() => {
        cy.visit('/');
    });

    it('clears all nodes and edges after clicking the "Clear" button', () => {
        // Click the "Clear" button
        cy.contains('button', 'Clear').click();
    
        // Verify no nodes exist in the ReactFlow pane
        cy.get('.react-flow__nodes')
            .find('.node')
            .should('not.exist');
    
        // Verify no edges exist in the ReactFlow pane
        cy.get('.react-flow__edges')
            .find('.react-flow__edge')
            .should('not.exist');
    });
});

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

describe('Drag and Drop Nodes', () => {
    beforeEach(() => {
      cy.visit('/');
    });
  
    it('drag and drop an input node onto the flowchart', () => {
      // Find the button containing "in" and drag it onto the ReactFlow pane
      cy.contains('button', 'in').dragAndDrop('.react-flow');
  
      // Verify the input node exists in the reactflow pane
      cy.get('.react-flow__nodes')
        .find('[data-testid=input-node]')
        .should('exist');
    });

    it('drag and drop an output node onto the flowchart', () => {
        // Find the button containing "out" and drag it onto the ReactFlow pane
        cy.contains('button', 'out').dragAndDrop('.react-flow');
    
        // Verify the input node exists in the reactflow pane
        cy.get('.react-flow__nodes')
            .find('[data-testid=output-node]')
            .should('exist');
    });

    it('drag and drop a binary node onto the flowchart', () => {
        // Find the button containing "+" and drag it onto the ReactFlow pane
        cy.contains('button', '+').dragAndDrop('.react-flow');
    
        // Verify the input node exists in the reactflow pane
        cy.get('.react-flow__nodes')
          .find('[data-testid=binary-node]')
          .should('exist');
    });
    
    it('drag and drop a unary node onto the flowchart', () => {
        // Find the button containing "sin" and drag it onto the ReactFlow pane
        cy.contains('button', 'sin').dragAndDrop('.react-flow');
    
        // Verify the input node exists in the reactflow pane
        cy.get('.react-flow__nodes')
        .find('[data-testid=unary-node]')
        .should('exist');
    });

    it('drag and drop a comparison node onto the flowchart', () => {
        // Find the button containing ">" and drag it onto the ReactFlow pane
        cy.contains('button', '>').dragAndDrop('.react-flow');
    
        // Verify the input node exists in the reactflow pane
        cy.get('.react-flow__nodes')
        .find('[data-testid=comparison-node]')
        .should('exist');
    });
});

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