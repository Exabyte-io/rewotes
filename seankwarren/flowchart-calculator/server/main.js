import { Meteor } from 'meteor/meteor';
import { FlowsCollection } from '/imports/api/flows';

Meteor.startup(async () => {
    // We publish the entire LFlowsinks collection to all clients.
    // In order to be fetched in real-time to the clients
    Meteor.publish("flows", function () {
        return FlowsCollection.find();
    });
});

Meteor.methods({

    saveFlow(data) {
        // Check if a flow with the same name already exists
        const existingFlow = FlowsCollection.findOne({ name: data.name });

        if (existingFlow) {
            // Update the existing flow with the new data (nodes and edges)
            FlowsCollection.update(existingFlow._id, {
                $set: {
                    nodes: data.nodes,
                    edges: data.edges,
                },
            });
        } else {
            // If the flow with the same name doesn't exist, insert the new flow
            FlowsCollection.insert(data);
        }
    },

    fetchFlows() {
        return FlowsCollection.find().fetch();
    },

    clearFlows() {
        FlowsCollection.remove({});
    },
});
