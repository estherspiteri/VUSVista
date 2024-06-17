import { useLocation } from "react-router-dom";
import PublicationViewPage from "../components/publication-view-page/publication-view-page";
import { publicationService } from "../services/publication/publication.service";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { IPublicationPreview } from "../models/publication-view.model";
import { IVUSSummary } from "../models/vus-summary.model";
import { convertPubDates } from "../helpers/date-helper";

const PublicationPhenotypeViewPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [publications, setPublications] =
    useState<IPublicationPreview[]>(undefined);
  const [variant, setVariant] = useState<IVUSSummary>(undefined);

  const loc = useLocation();
  const variantId = loc.pathname.split("/publication-phenotype-view/")[1];
  const searchParams = new URLSearchParams(loc.search);
  const rsid = searchParams.get("rsid");
  const optionalText = searchParams.get("phenotype");

  useEffect(() => {
    if (isLoading) {
      publicationService
        ?.getPublicationsByRsidAndWithOptionalText({
          variantId: variantId,
          rsid: rsid,
          optionalText: optionalText,
        })
        .then((res) => {
          if (res.isSuccess) {
            setPublications(convertPubDates(res.publications));
            setVariant(res.variant);
            setIsLoading(false);
          } else {
            //TODO: handle error
          }
        });
    }
  }, [isLoading, rsid, optionalText]);

  if (isLoading) {
    return <Loader />;
  } else {
    return (
      <PublicationViewPage
        isPhenotypePublicationPage={true}
        variantId={variantId}
        publications={publications}
        variant={variant}
        phenotype={optionalText}
      />
    );
  }
};

export default PublicationPhenotypeViewPageWrapper;
