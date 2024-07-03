import React, { createContext, useEffect, useState } from "react";
import { vusService } from "./services/vus/vus.service";
import { IStatus } from "./services/vus/vus.dto";

// Initial context value
const initialContext = {
  taskIds: [],
  setTaskIds: (updatedIds: number[]) => {},
  completedTasks: [],
  setCompletedTasks: (updatedStatuses: IStatus[]) => {},
  isUserLoggedIn: false,
  setIsUserLoggedIn: (isLoggedIn: boolean) => {},
};

export const AppContext = createContext(initialContext);

export const AppProvider = ({ children }) => {
  const [taskIds, setTaskIds] = useState<number[]>(initialContext.taskIds);
  const [completedTasks, setCompletedTasks] = useState<IStatus[]>(
    initialContext.completedTasks
  );
  const [isUserLoggedIn, setIsUserLoggedIn] = useState(
    initialContext.isUserLoggedIn
  );
  console.log(taskIds);
  // poll to check if any file upload tasks have completed
  useEffect(() => {
    const pollTaskStatus = () => {
      const interval = setInterval(async () => {
        if (taskIds.length === 0) return;

        vusService
          .checkFileUploadStatuses({
            taskIds: taskIds.map((i) => i.toString()),
          })
          .then((res) => {
            const completedStatuses = res.statuses?.filter(
              (s) => s.isSuccess !== null
            );

            setCompletedTasks(completedTasks.concat(completedStatuses));

            let updatedTaskIds = taskIds.filter(
              (id) => !completedStatuses.map((s) => s.taskId).includes(id)
            );

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
    <AppContext.Provider
      value={{
        taskIds,
        setTaskIds,
        isUserLoggedIn,
        setIsUserLoggedIn,
        completedTasks,
        setCompletedTasks,
      }}
    >
      {children}
    </AppContext.Provider>
  );
};
