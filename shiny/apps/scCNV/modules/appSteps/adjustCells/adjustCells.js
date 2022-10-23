cellKeepReject = function(prefix, cell_id, override){
    Shiny.setInputValue(
        prefix + 'cellKeepReject',
        {cell_id: cell_id, override: override},
        {priority: "event"}
    );
};
cellSetModalCN = function(prefix, cell_id, inputId){
    Shiny.setInputValue(
        prefix + 'cellSetModalCN',
        {cell_id: cell_id, value: $("#" + inputId).val()},
        {priority: "event"}
    );
};
