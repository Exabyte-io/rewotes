describe('Clear Nodes and Edges', () => {
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