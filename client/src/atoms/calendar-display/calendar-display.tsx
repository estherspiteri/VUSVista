import React from "react";
import "./calendar-display.scss";
import Calendar from "react-calendar";

type CalendarDisplayProps = {
  markedDates?: { date: string; update: boolean }[];
};

const CalendarDisplay: React.FunctionComponent<CalendarDisplayProps> = (
  props: CalendarDisplayProps
) => {
  return (
    <div className="calendar">
      <Calendar
        tileClassName={({ date, view }) => {
          const formattedDate = `${date.getFullYear()}/${(date.getMonth() + 1)
            .toString()
            .padStart(2, "0")}/${date.getDate().toString().padStart(2, "0")}`;

          const markedDate = props.markedDates.find((x) => {
            return x.date === formattedDate;
          });
          if (view !== "month") return null;
          if (markedDate) {
            if (markedDate.update) {
              return "highlight-update";
            }
            return "highlight";
          }
        }}
      />
    </div>
  );
};

export default CalendarDisplay;
