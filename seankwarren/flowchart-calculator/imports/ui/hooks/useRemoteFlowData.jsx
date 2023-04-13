import React, { useCallback, useState } from 'react';
import useInputField from './useInputField';

const useRemoteFlowData = ({ nodes, edges }) => {
    // state for fetching saved flows from db
    const [fetchedFlows, setFetchedFlows] = useState([]);
    const { flowName, updateFlowName } = useInputField();

    function updateFetchedFlows(flows) {
        setFetchedFlows(flows);
    }

    const saveFlow = useCallback(() => {
        // Check if a flow with the same name already exists
        const flowWithNameExists = fetchedFlows.some(
            (flow) => flow.name === flowName
        );

        if (flowWithNameExists) {
            alert(
                'A flow with this name already exists. Please choose a different name.'
            );
        } else {
            console.log(nodes, edges);
            Meteor.call(
                'saveFlow',
                { nodes, edges, name: flowName },
                (error) => {
                    if (error) {
                        console.log(error.reason);
                    } else {
                        // updateFlowName(''); // Clear the flow name input field after a successful save
                        fetchFlows(); // Fetch the flows again to update the dropdown
                    }
                }
            );
        }
    }, [nodes, edges]);

    const fetchFlows = () => {
        Meteor.call('fetchFlows', (error, result) => {
            if (error) {
                console.error('Error fetching flows:', error);
            } else {
                updateFetchedFlows(result);
            }
        });
    };

    const clearFlows = () => {
        Meteor.call('clearFlows', (error) => {
            if (error) {
                console.error('Error clearing flows:', error);
            } else {
                console.log('Flows collection cleared');
            }
        });
    };

    return {
        fetchedFlows,
        saveFlow,
        fetchFlows,
        clearFlows,
        flowName,
        updateFlowName,
    };
};

export default useRemoteFlowData;
