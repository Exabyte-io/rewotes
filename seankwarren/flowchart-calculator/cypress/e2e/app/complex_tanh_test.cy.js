/// <reference types="cypress" />

describe('Complex Flow test (tanh)', () => {
    beforeEach(() => {
        // Visit the app's home page before each test
        cy.visit('/');
    });
  
    it('tests a complex flowchart', () => {

        function addNode(nodeTypeLabel, nodeIdArray) {
            return cy.contains('button', nodeTypeLabel).dragAndDrop('.react-flow').then(() => {
                const nodeTypes = {
                    'in': 'input',
                    'out': 'output',
                    '>': 'comparison',
                    '+': 'binary',
                    'sin': 'unary',
                }
                console.log(`${nodeTypes[nodeTypeLabel]}-node`);
                return cy.get(`[data-testid=${nodeTypes[nodeTypeLabel]}-node]`).then(($nodes) => {
                    $nodes.each((_, el) => {
                        const currentNodeId = el.getAttribute('data-nodeid');
                        if (!nodeIdArray.includes(currentNodeId)) {
                            nodeIdArray.push(currentNodeId);
                            console.log(currentNodeId, nodeIdArray)
                            return false;
                        }
                    });
                });
            });
        }

        const inputNodeIds = [];
        const unaryNodeIds = [];
        const binaryNodeIds = [];
        const outputNodeIds = [];

        // Click the clear button
        cy.get('.clear').click();
    
        // Place new nodes on the flowchart and update nodeId arrays
        addNode('in', inputNodeIds)
            .then(() => addNode('in', inputNodeIds))
            .then(() => addNode('in', inputNodeIds))
            .then(() => addNode('sin', unaryNodeIds))
            .then(() => addNode('+', binaryNodeIds))
            .then(() => addNode('+', binaryNodeIds))
            .then(() => addNode('+', binaryNodeIds))
            .then(() => addNode('+', binaryNodeIds))
            .then(() => addNode('out', outputNodeIds))
            .then(() => {
                cy.get(`[data-nodeid=${unaryNodeIds[0]}] select`)
                    .select('e^x', { force: true });

                cy.get(`[data-nodeid=${binaryNodeIds[0]}] select`)
                    .select('^', { force: true });
                
                cy.get(`[data-nodeid=${binaryNodeIds[1]}] select`)
                    .select('-', { force: true });
                
                cy.get(`[data-nodeid=${binaryNodeIds[2]}] select`)
                    .select('+', { force: true });

                cy.get(`[data-nodeid=${binaryNodeIds[3]}] select`)
                    .select('/', { force: true });

                // Connect the nodes
                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${inputNodeIds[0]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${unaryNodeIds[0]}] .handle.target`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${unaryNodeIds[0]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[0]}] .handle.target.top`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${inputNodeIds[1]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[0]}] .handle.target.bottom`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${inputNodeIds[2]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[1]}] .handle.target.bottom`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${inputNodeIds[2]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[2]}] .handle.target.bottom`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[0]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[1]}] .handle.target.top`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[0]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[2]}] .handle.target.top`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[1]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[3]}] .handle.target.top`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[2]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[3]}] .handle.target.bottom`
                );

                cy.connectHandles(
                    `.react-flow__nodes [data-nodeid=${binaryNodeIds[3]}] .handle.source`,
                    `.react-flow__nodes [data-nodeid=${outputNodeIds[0]}] .handle.target`
                );


                // Set input values and check the output value
                cy.get(`[data-nodeid=${inputNodeIds[0]}] input`)
                    .clear({ force: true })
                    .type('3', { force: true });

                cy.get(`[data-nodeid=${inputNodeIds[1]}] input`)
                    .clear({ force: true })
                    .type('2', { force: true });
                
                cy.get(`[data-nodeid=${inputNodeIds[2]}] input`)
                    .clear({ force: true })
                    .type('1', { force: true });

                cy.get(`[data-nodeid=${outputNodeIds[0]}]`)
                    .should('contain.text', '0.99505475368');
            });
    
    });
});