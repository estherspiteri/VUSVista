import React, { useState } from "react";
import "./calendar-display.scss";
import Calendar from "react-calendar";

type CalendarDisplayProps = {
  markedDates?: { date: string; update: boolean }[];
};

const CalendarDisplay: React.FunctionComponent<CalendarDisplayProps> = (
  props: CalendarDisplayProps
) => {
  const [activeStartDate, setActiveStartDate] = useState(new Date());

  return (
    <div className="calendar">
      <Calendar
        tileClassName={({ date, view }) => {
          if (view !== "month") return null;
          const startOfMonth = new Date(
            activeStartDate.getFullYear(),
            activeStartDate.getMonth(),
            1
          );
          const endOfMonth = new Date(
            activeStartDate.getFullYear(),
            activeStartDate.getMonth() + 1,
            0
          );

          //hide dates from previous/next month
          if (date < startOfMonth || date > endOfMonth) {
            return "hide";
          }

          const formattedDate = `${date.getFullYear()}/${(date.getMonth() + 1)
            .toString()
            .padStart(2, "0")}/${date.getDate().toString().padStart(2, "0")}`;

          const markedDate = props.markedDates.find((x) => {
            return x.date === formattedDate;
          });

          //highlight marked dates & make bold those dates with an update
          if (markedDate) {
            if (markedDate.update) {
              return "highlight-update";
            }
            return "highlight";
          }
        }}
        tileContent={tileContent}
        onActiveStartDateChange={(res) =>
          setActiveStartDate(res.activeStartDate)
        }
      />
    </div>
  );

  function tileContent({ date, view }) {
    if (view === "month") {
    }
    return null;
  }
};

export default CalendarDisplay;
