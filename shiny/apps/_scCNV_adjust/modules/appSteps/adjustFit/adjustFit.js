setAdjustFitCell = function(prefix, cell_id){
    Shiny.setInputValue(
        prefix + 'setAdjustFitCell',
        cell_id,
        {priority: "event"}
    );
};
