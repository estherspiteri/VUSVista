import React, { createContext, useEffect, useState } from "react";
import { vusService } from "./services/vus/vus.service";

// Initial context value
const initialContext = {
  taskIds: [],
  setTaskIds: (updatedIds: number[]) => {},
};

export const AppContext = createContext(initialContext);

export const AppProvider = ({ children }) => {
  const [taskIds, setTaskIds] = useState(initialContext.taskIds);
  const [pollingInterval, setPollingInterval] = useState(null);

  useEffect(() => {
    if (taskIds.length === 0) return;

    const pollTaskStatus = () => {
      const interval = setInterval(async () => {
        vusService.checkFileUploadStatuses({ taskIds: taskIds }).then((res) => {
          let updatedTaskIds = taskIds;

          res.statuses?.forEach((status) => {
            if (status.isSuccess !== null) {
              updatedTaskIds = updatedTaskIds.filter(
                (id) => id !== status.taskId
              );

              if (status.isSuccess) {
                console.log("SUCCESS");
              } else {
                console.log("FAILLLL");
              }
            }
          });
          setTaskIds(updatedTaskIds);

          if (updatedTaskIds.length === 0) {
            clearInterval(interval);
          }
        });
      }, 30000);

      return () => clearInterval(interval);
    };

    const intervalCleanup = pollTaskStatus();

    return () => {
      if (intervalCleanup) {
        intervalCleanup();
      }
    };
  }, [taskIds]);

  return (
    <AppContext.Provider value={{ taskIds, setTaskIds }}>
      {children}
    </AppContext.Provider>
  );
};
