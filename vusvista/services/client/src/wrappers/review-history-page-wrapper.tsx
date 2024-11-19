import { useLocation } from "react-router-dom";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { reviewService } from "../services/review/review.service";
import { IVUSSummary } from "../models/vus-summary.model";
import ReviewHistoryPage from "../components/review-history-page/review-history-page";
import { IClassificationReview } from "../models/classification-review.model";
import { convertReviewDates } from "../helpers/date-helper";

const ReviewPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [variantSummary, setVariantSummary] = useState<IVUSSummary>(undefined);
  const [reviews, setReviews] = useState<IClassificationReview[]>(undefined);

  const loc = useLocation();
  const vusId = loc.pathname.split("/review-history/")[1];

  useEffect(() => {
    if (isLoading) {
      reviewService
        .getAllClassificationReviews({ vusId: parseInt(vusId) })
        .then((res) => {
          setVariantSummary(res.variantSummary);
          setReviews(convertReviewDates(res.reviews));
          setIsLoading(false);
        });
    }
  }, [isLoading, vusId]);

  if (isLoading) {
    return <Loader />;
  } else {
    return (
      <ReviewHistoryPage
        variantSummary={variantSummary}
        reviews={reviews}
      />
    );
  }
};

export default ReviewPageWrapper;
