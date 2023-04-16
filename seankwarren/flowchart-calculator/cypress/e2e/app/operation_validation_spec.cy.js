/// <reference types="cypress" />

describe('Binary Operators', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('tests the binary operators', () => {
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
            .type('3', { force: true }); // Set the new value to '42'
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '6');
        
        cy.get('[data-testid=binary-node] select')
            .select('-', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '0');

        cy.get('[data-testid=binary-node] select')
            .select('*', { force: true });

        cy.get('[data-testid=output-node]')
            .should('contain.text', '9');

        cy.get('[data-testid=binary-node] select')
            .select('/', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '1');
        
        cy.get('[data-testid=binary-node] select')
            .select('^', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '27');
    });
});

describe('Unary Operators', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('tests the unary operators', () => {
        // Click the clear button
        cy.get('.clear').click();
    
        // Place a new node on the flowchart
        cy.contains('button', 'in').dragAndDrop('.react-flow');
        cy.contains('button', 'sin').dragAndDrop('.react-flow');
        cy.contains('button', 'out').dragAndDrop('.react-flow');
    
        // Connect the nodes
        cy.connectHandles(
            '.react-flow__nodes [data-testid=input-node] .handle.source',
            '.react-flow__nodes [data-testid=unary-node] .handle.target'
        );

        cy.connectHandles(
            '.react-flow__nodes [data-testid=unary-node] .handle.source',
            '.react-flow__nodes [data-testid=output-node] .handle.target'
        );
    
        cy.get('[data-testid=input-node] input')
            .clear({ force: true }) // Clear the current value of the input element
            .type('3.14159', { force: true }); // Set the new value to '42'
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '0.00000265358979335273');
        
        cy.get('[data-testid=unary-node] select')
            .select('cos', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '-0.9999999999964793')
        
        cy.get('[data-testid=unary-node] select')
            .select('tan', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '-0.000002653589793362073');
        
        cy.get('[data-testid=unary-node] select')
            .select('e^x', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', '23.14063122695496')
    });
});

describe('Comparison Operators', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('tests the comparison operators', () => {
        // Click the clear button
        cy.get('.clear').click();
    
        // Place a new node on the flowchart
        cy.contains('button', 'in').dragAndDrop('.react-flow');
        cy.contains('button', '>').dragAndDrop('.react-flow');
        cy.contains('button', 'out').dragAndDrop('.react-flow');
    
        // Connect the nodes
        cy.connectHandles(
            '.react-flow__nodes [data-testid=input-node] .handle.source',
            '.react-flow__nodes [data-testid=comparison-node] .handle.target.top'
        );

        cy.connectHandles(
            '.react-flow__nodes [data-testid=input-node] .handle.source',
            '.react-flow__nodes [data-testid=comparison-node] .handle.target.bottom'
        );

        cy.connectHandles(
            '.react-flow__nodes [data-testid=comparison-node] .handle.source',
            '.react-flow__nodes [data-testid=output-node] .handle.target'
        );
    
        cy.get('[data-testid=input-node] input')
            .clear({ force: true }) // Clear the current value of the input element
            .type('3.14159', { force: true }); // Set the new value to '42'
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', 'false');
        
        cy.get('[data-testid=comparison-node] select')
            .select('<', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', 'false');
        
        cy.get('[data-testid=comparison-node] select')
            .select('>=', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', 'true');
        
        cy.get('[data-testid=comparison-node] select')
            .select('<=', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', 'true');
        
        cy.get('[data-testid=comparison-node] select')
            .select('!=', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', 'false');
        
        cy.get('[data-testid=comparison-node] select')
            .select('==', { force: true });
        
        cy.get('[data-testid=output-node]')
            .should('contain.text', 'true')
    });
});