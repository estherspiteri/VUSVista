import { ILoadAllSamplesResponse } from "./sample.dto";

export class SampleService {
  async loadAllSamples(): Promise<ILoadAllSamplesResponse> {
    const result: ILoadAllSamplesResponse = await fetch(`/sample/view`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    let updatedSampleList = result?.sampleList?.map((s) => {
      return { ...s, dateOfFileUpload: new Date(s.dateOfFileUpload) };
    });

    return { isSuccess: result.isSuccess, sampleList: updatedSampleList };
  }
}

export const samplesService = new SampleService();
