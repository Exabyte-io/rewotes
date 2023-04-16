describe('Drag and Drop Input Node', () => {
    beforeEach(() => {
      cy.visit('/');
    });
  
    it('drags an input node onto the reactflow pane and checks if it exists', () => {
      // Find the button containing "in" and drag it onto the ReactFlow pane
      cy.contains('button', 'in').dragAndDrop('.react-flow');
  
      // Verify the input node exists in the reactflow pane
      cy.get('.react-flow__nodes')
        .find('[data-testid=input-node]')
        .should('exist');
    });
});

describe('Drag and Drop Output Node', () => {
    beforeEach(() => {
      cy.visit('/');
    });
  
    it('drags an output node onto the reactflow pane and checks if it exists', () => {
      // Find the button containing "in" and drag it onto the ReactFlow pane
      cy.contains('button', 'out').dragAndDrop('.react-flow');
  
      // Verify the input node exists in the reactflow pane
      cy.get('.react-flow__nodes')
        .find('[data-testid=output-node]')
        .should('exist');
    });
});

describe('Drag and Drop Binary Node', () => {
    beforeEach(() => {
      cy.visit('/');
    });
  
    it('drags a binary node onto the reactflow pane and checks if it exists', () => {
      // Find the button containing "in" and drag it onto the ReactFlow pane
      cy.contains('button', '+').dragAndDrop('.react-flow');
  
      // Verify the input node exists in the reactflow pane
      cy.get('.react-flow__nodes')
        .find('[data-testid=binary-node]')
        .should('exist');
    });
});

describe('Drag and Drop Unary Node', () => {
    beforeEach(() => {
      cy.visit('/');
    });
  
    it('drags a unary node onto the reactflow pane and checks if it exists', () => {
      // Find the button containing "in" and drag it onto the ReactFlow pane
      cy.contains('button', 'sin').dragAndDrop('.react-flow');
  
      // Verify the input node exists in the reactflow pane
      cy.get('.react-flow__nodes')
        .find('[data-testid=unary-node]')
        .should('exist');
    });
});

describe('Drag and Drop Comparison Node', () => {
    beforeEach(() => {
      cy.visit('/');
    });
  
    it('drags a comparison node onto the reactflow pane and checks if it exists', () => {
      // Find the button containing "in" and drag it onto the ReactFlow pane
      cy.contains('button', '>').dragAndDrop('.react-flow');
  
      // Verify the input node exists in the reactflow pane
      cy.get('.react-flow__nodes')
        .find('[data-testid=comparison-node]')
        .should('exist');
    });
});
